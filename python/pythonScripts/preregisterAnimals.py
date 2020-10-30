import SimpleITK as sitk
import os
import argparse
import numpy as np
from glob import glob
import tifffile
import javabridge
import bioformats
import io, sys
import rotate, rotate2
from sklearn.decomposition import PCA
from tqdm import tqdm

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def resample(image, transform, factor=1.0):
    # Output image Origin, Spacing, Size, Direction are taken from the reference
    # image in this call to Resample
    #reference_image=image
    reference_image = np.zeros((int(factor*image.GetSize()[2]), int(factor*image.GetSize()[1]), int(factor*image.GetSize()[0]))) 
    reference_image = sitk.GetImageFromArray(reference_image)
    interpolator = sitk.sitkCosineWindowedSinc
    default_value = 0.0 # TODO
    return sitk.Resample(image, reference_image, transform,
                         interpolator, default_value)

def affine_rotate_3d(grid, transform,  rotation=None, affine_center=None, offset=None, factor=1.0):
    new_transform = sitk.AffineTransform(transform)
    matrix = np.array(transform.GetMatrix()).reshape((3, 3)) 
    if rotation is None:
        print("No rotation observed")
        new_matrix = matrix
        new_transform.SetMatrix(new_matrix.ravel())
    else:
        new_matrix = np.dot(rotation, matrix)
        new_transform.SetMatrix(new_matrix.ravel())
        new_transform.SetCenter((affine_center[2], affine_center[1], affine_center[0]))
    
    if offset is None:
        pass
    else:
        new_transform.SetTranslation((offset[2], offset[1], offset[0] )) # xyz
    #new_transform.SetTranslation((0, 0, 0))
    resampled = resample(grid, new_transform, factor=factor)
    print(new_matrix)
    return resampled, new_transform


def getIsotropicVersions(image, pixel_size):
    image=np.transpose(image, (0,2,3,1)) # channel is last dim
    image_resampled=[]
    for ch in range(image.shape[-1]):
        image_ch=image[..., ch]
        image_ch_sitk=sitk.GetImageFromArray(image_ch)
        image_ch_sitk.SetSpacing(pixel_size)
        image_ch_sitk_resampled=make_isotropic(image_ch_sitk)
        image_ch_resampled = sitk.GetArrayFromImage(image_ch_sitk_resampled)
        image_resampled.append(image_ch_resampled)
    return np.asarray(image_resampled)

def getPixelSize(imagePath):
    javabridge.start_vm(class_path=bioformats.JARS)
    path_to_data = imagePath
    xml_string = bioformats.get_omexml_metadata(path_to_data)
    #ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
    pos_z=xml_string.find("PhysicalSizeZ=")
    pos_y=xml_string.find("PhysicalSizeY=")
    pos_x=xml_string.find("PhysicalSizeX=")
    return np.array([float(xml_string[pos_x+15:pos_x+20]), float(xml_string[pos_y+15:pos_y+20]), float(xml_string[pos_z+15:pos_z+20])]) 

def make_isotropic(image, interpolator = sitk.sitkLinear):
    '''
    Resample an image to isotropic pixels (using smallest spacing from original) and save to file. Many file formats 
    (jpg, png,...) expect the pixels to be isotropic. By default the function uses a linear interpolator. For
    label images one should use the sitkNearestNeighbor interpolator so as not to introduce non-existant labels.
    '''
    original_spacing = image.GetSpacing()
    # Image is already isotropic, just return a copy.
    if all(spc == original_spacing[0] for spc in original_spacing):
        return sitk.Image(image)
    # Make image isotropic via resampling.
    original_size = image.GetSize()
    min_spacing = min(original_spacing)
    new_spacing = [min_spacing]*image.GetDimension()
    new_size = [int(round(osz*ospc/min_spacing)) for osz,ospc in zip(original_size, original_spacing)]
    return sitk.Resample(image, new_size, sitk.Transform(), interpolator,
                         image.GetOrigin(), new_spacing, image.GetDirection(), 0,
                         image.GetPixelID())


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ch_prereg", default=0, type=int,
            help="channel on which pca/detections should be done ")
    parser.add_argument("--input_directory_name", default="/projects/fileserver/projects/Project_Johannes/Imaging/Neuropeptides_Ec_Pc_Ml/Neuropeptides/prereg/prereg_test/", type=str, help="input directory name")
    parser.add_argument("--output_directory_name", default="/projects/fileserver/projects/Project_Johannes/Imaging/Neuropeptides_Ec_Pc_Ml/Neuropeptides/prereg/prereg_test_output/", type=str, help="output directory name");
    parser.add_argument("--dsFactor", default=5, type=int, help="down sampling factor");
    
    args = parser.parse_args();
    return args;

def main():
    args = parse_args()
    dsFactor=args.dsFactor
    ch_prereg=args.ch_prereg
    filenames = sorted(glob(args.input_directory_name+'*.tif'))
    
    for filename in tqdm(filenames[:1]):
        image = tifffile.imread(filename)
        print("Image from ", os.path.basename(filename), " is of shape",image.shape)         
        text_trap=io.StringIO()
        sys.stdout=text_trap
        image_pixel_size=getPixelSize(filename)
        sys.stdout=sys.__stdout__
        print("Image from ", os.path.basename(filename), " has pixel sizes of ",image_pixel_size)         
        image = getIsotropicVersions(image, image_pixel_size)
        print("Isotropic Image from ", os.path.basename(filename), " is of shape",image.shape)         
        image =np.transpose(image, (1, 0, 2, 3))
        image_channel = image[:, ch_prereg, ...][::dsFactor, ::dsFactor, ::dsFactor].astype('float32')
        d1, h1, w1= image_channel.shape
        image_channel_expanded=np.zeros((int(1.2*image_channel.shape[0]), int(1.2*image_channel.shape[1]), int(1.2*image_channel.shape[2])))
        d_new, h_new, w_new =  image_channel_expanded.shape
        image_channel_expanded[d_new//2-d1//2:d_new//2-d1//2+d1, h_new//2-h1//2:h_new//2-h1//2+h1, w_new//2-w1//2 :w_new//2-w1//2+w1] = image_channel
        image_channel=image_channel_expanded
        detections=rotate.findBlobs(image_channel, scales=range(int(5/dsFactor), int(10/dsFactor)), threshold = 10)
        pca=PCA(n_components=3)
        pca.fit(np.flip(detections[:, 1:], 1)) # TODO
        V=pca.components_
        rigid_transform = np.linalg.inv(V)
        rotation_matrix = rigid_transform
        image_channel_sitk=sitk.GetImageFromArray(image_channel)
        centerEmbryo=np.mean(detections[:,1:], 0)
        dimension = 3 
        affine = sitk.AffineTransform(dimension)
        image_channel_sitk_warped, affine = affine_rotate_3d(image_channel_sitk, affine, np.linalg.inv(rotation_matrix), centerEmbryo, offset=None, factor=1.2) #TODO
        image_channel_sitk_warped_array = sitk.GetArrayFromImage(image_channel_sitk_warped)
        detectionszyx = detections[:, 1:] # zyx
        detectionszyx_copy=np.zeros_like(detectionszyx)
        detectionszyx_copy[:, 0]=detectionszyx[:, 2] # x--> z
        detectionszyx_copy[:, 1]=detectionszyx[:, 1]
        detectionszyx_copy[:, 2]=detectionszyx[:, 0] # z--> x
        detectionszyx_copy = detectionszyx_copy.transpose()
        detectionszyx_copy=detectionszyx_copy-np.flip(centerEmbryo[:, np.newaxis]) # zyx
        transformedDetections=np.matmul(rotation_matrix, detectionszyx_copy) + np.flip(centerEmbryo[:, np.newaxis]) 
        transformedDetections= transformedDetections.transpose() # --> xyz
        centerEmbryoTransformed=np.flip(np.mean(transformedDetections, 0))
        expandedImageCenter=np.array([image_channel_sitk_warped_array.shape[0]//2, image_channel_sitk_warped_array.shape[1]//2, image_channel_sitk_warped_array.shape[2]//2])
        correctOffset=expandedImageCenter-centerEmbryoTransformed # --> zyx
        affine = sitk.AffineTransform(dimension)
        image_channel_sitk_warped_2, affine = affine_rotate_3d(image_channel_sitk_warped, affine, rotation=None, offset=-correctOffset, factor=1.0)
        image_channel_sitk_warped_array_2 = sitk.GetArrayFromImage(image_channel_sitk_warped_2)
        transformedDetections=transformedDetections+np.flip(correctOffset)
        brightDetections=rotate.findBlobs(image_channel_sitk_warped_array_2, scales=range(int(5/dsFactor), int(10/dsFactor)), threshold = 10)
        intensityLeft=0
        intensityRight=0
        countRight=0
        countLeft=0
        for detection in brightDetections:
            if detection[3] > expandedImageCenter[2]:
                intensityRight+=image_channel_sitk_warped_array_2[detection[1],detection[2], detection[3]]
                countRight+=1
            else:
                intensityLeft+=image_channel_sitk_warped_array_2[detection[1],detection[2], detection[3]]
                countLeft+=1
        if(intensityLeft > intensityRight):
            print("head is currently on the right: head should be flipped!")
            image_channel_sitk_warped_array_3=np.rot90(image_channel_sitk_warped_array_2, 2, axes=(1, 2)) # TODO
        else:
            print("head is on the left! no need for flipping!")
            image_channel_sitk_warped_array_3=image_channel_sitk_warped_array_2


        #darkDetections=rotate2.findBlobs(image_channel_sitk_warped_array_2, scales=range(int(20/dsFactor), int(40/dsFactor)), threshold = 4)
        #indices1, =np.where(darkDetections[:, 3] > expandedImageCenter[2])
        #indices2, =np.where(darkDetections[:, 3] <= expandedImageCenter[2])
        #print("Indices", len(indices1), len(indices2))
        #if(len(indices1) < len(indices2)):
        #    print("head is currently on the right: head should be flipped!")
        #    image_channel_sitk_warped_array_3=np.rot90(image_channel_sitk_warped_array_2, 2, axes=(1, 2)) # TODO
        #else:
        #    print("head is on the left! no need for flipping!")
        #    image_channel_sitk_warped_array_3=image_channel_sitk_warped_array_2
        tifffile.imsave(args.output_directory_name+os.path.basename(filename)[:-4]+"_prereg.tif", image_channel_sitk_warped_array_3.astype(np.float32))
        #tifffile.imsave(args.output_directory_name+os.path.basename(filename)[:-4]+"_prereg.tif", np.max(image_channel_sitk_warped_array_3, 0).astype(np.float32))
	
if __name__ == "__main__":
    main();
