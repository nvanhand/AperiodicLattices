# AperiodicLattices
## File Types
For the project, there are 2 primary data types - images and STLs. 
### Images
Images represent lattices with pixels. Based on the file type and how it was saved, each pixel can represent a range of values and may have as many as 3 values associated with each pixel, as in the case of RGB. In general, we can consider an image an array of values, with the depth representing dimensions of color. For our purposes, we typically want to simplify images to "lattice" and "not lattice" requiring the binarization of each image so the arrays are 1D. Images do not exist in "real" measurements - we almost must translate between the image dimensions and some sort of scale. In some cases, if we don't trust our lens or the way we took the picture, we may need to also perform correction for lens distortion. We'll talk about this another time. 

In images, we typically measure features in number of pixels and use a multiplicative conversion to get real values. For the Keyence images (PNGs), I just manually measured the scalebars - all were about 213 px per cm (+/- 1 px). For TIFs, i used a slightly more complex process as the resolution was higher - see tif_to_scale.py. 

### Standard Triangle Language (STL)
STLs are saved in real values. Rather than an array of values, STLs are closer to lists of vertices and edges. These coordinates all exist in space which relate them. As a result, we do not 

## Basic Functions
