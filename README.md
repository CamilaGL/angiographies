# angiographies

Processing angiographies.

# Conda environments

## angiographies

### Created with:
conda create --name angiographies<br>
conda install python=3.6 scikit-image simpleitk vmtk -c simpleitk -c vmtk<br>
conda install llvm=3.3<br>
pip install -e . #to install project <br>

## edt

### Created with:
edt (euclidean distance transform)<br>
conda create –name edt<br>
conda install edt -c conda-forge python=3.9<br>
conda install simpleitk -c simpleitk<br>
conda install vtk<br>
conda install vmtk -c conda-forge<br>
pip install -e . #to install project <br>

## edt2

### Created with:<br>
conda create --name edt2 --clone edt<br>
conda install scikit-image -c conda-forge<br>


# Usage

## Thinning / skeletonisation

## vmtknetworkextraction

## Nidus extraction


## Authors
Camila García

