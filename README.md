# Practices 2.10, 2.11, 2.12, and 2.13

## Practice 2.10 

### Thermal Equation with bilinear quadrilateral elemets 

<!--
### FEM2D: Linear Quadrilateral FiniteElements
-->
## Practice 2.11

### Thermal Equation with bilinear quadrilateral elemets and convection B.C.

<!--
### FEM2D: Thermal convection problems using Quadrilateral Elements
-->

## Practice 2.12 

### Convection BC example: triangular version

<!--
### FEM2D: Thermal convection problems using Triangular Elements

#### Engine cooling fin
-->


## Practice 2.13 

### The same problem using a triangular and quadrilateral elements and convection BC

<!--
FEM2D: Thermal convection problems using Quadrilateral elements

#### Engine cooling fin
-->

## Materials
The guides that explain these practices are available at prof. Antonio Susín's [Numerical Factory: FEM -- Practices](https://numfactory.upc.edu/numfactory/subjects/?subject=FEM&version=etseib&lang=en&section=2&subsection=P)

* [Guide of practice 10](https://numfactory.upc.edu/numfactory/subjects/FEM/FEM/practices/P2.9/html/P2_9_FEM2DThermalQuad_en.html)
* [Guide of practice 11](https://numfactory.upc.edu/numfactory/subjects/FEM/FEM/practices/P2.10/html/P2_10_FEM2DThermalQuadConvection_en.html)
* [Guide of practice 12](https://numfactory.upc.edu/numfactory/subjects/FEM/FEM/practices/P2.11a/html/P2_11a_FEM2DAletaThermalTriangConvection_en.html)
* [Guide of practice 13](https://numfactory.upc.edu/numfactory/subjects/FEM/FEM/practices/P2.11b/html/P2_11b_FEM2DAletaThermalQuadConvection_en.html)

The required extra and mesh files are in two zipped folders that can be downloaded from the same web page (click the links you will find at the bottom). We also give the links here:

* [addtionalFiles.zip](https://numfactory.upc.edu/web/FiniteElements/files/additionalFiles.zip)
* [meshFilesAll.zip](https://numfactory.upc.edu/numfactory/subjects/FEM/Extra/meshFilesAll.zip)

## Extra and plot files required

* `linearTriangElment.m`
* `applyConvTriang.m`
* `baryCoord.m`  
* `gaussValues1D.m`
* `gaussValues2DQuad.m`
* `bilinearQuadElement.m`
* `applyConvQuad.m`
* `baryCoordQuad.m`
* `plotElements.m`
* `plotElementsOld.m`
* `plotCountourSolution.m`

Uncompress the zipped folder: [addtionalFiles.zip](https://numfactory.upc.edu/web/FiniteElements/files/additionalFiles.zip) and place these files in in the same directory that the scripts which call the corresponding functions. 

Additionally, the following mesh files are required as well:

* `mesh2x2Quad.m`
* `mesh4x4Quad.m`
* `mesh8x8Quad.m`
* `meshmeshPlacaForatQuad.m`
* `meshAleta2DQuad.m`
* `meshAleta2DTriang.m`

Uncompress the zipped folder: [addtionalFiles.zip](https://numfactory.upc.edu/numfactory/subjects/FEM/Extra/meshFilesAll.zip) and place these files in in the same directory that the scripts which load the corresponding mesh files.  

#### Disclaimer

This stuff is provided 'as is'. Please, check it (if you find it useful), but try to write the programs on your own. In any case, using the content of  this repository is under your responsibility. 

If you find any mistakes (or have any suggestions), please report them to 

juan.ramon.pacha@upc.edu 

Many thanks,

J.R.
