import SimpleITK as sitk
import argparse
import os
from collections import OrderedDict
import json
import hashlib
from datetime import datetime



# def getDICOMMetadataVTK(dicomdir):
#     #read DICOM metadata with VTK
#     #tags: https://www.dicomlibrary.com/dicom/dicom-tags/
#     dparser = vtk.vtkDICOMParser()
#     dparser.SetDirectoryName(dicomdir)
#     dparser.Update()
    
#     metaData = dparser.GetMetaData()

#     if metaData.Has(vtk.vtkDICOMTag(0x0010, 0x0020)):
#         print (metaData.GetAttributeValue(vtk.vtkDICOMTag(0x0010, 0x0020)))


metaIDs = OrderedDict()
metaIDs["patientID"] = '0010|0020'
metaIDs["studyDate"] = '0008|0020'
metaIDs["studyID"] = '0020|0010'

metaKeys = OrderedDict()
metaKeys['0010|0020'] = "patientID"
metaKeys['0020|0010'] = "studyID"
metaKeys['0008|0020'] = "study date"
metaKeys['0020|000e'] = "seriesID"
#metaKeys['0008|0080'] = "institution name"
metaKeys['0010|0030'] = "patient birth date"
metaKeys['0010|0040'] = "patient assigned sex"
metaKeys['0010|1010'] = "patient age"
#metaKeys['0028|0006'] = "image dimensions"
#metaKeys['0028|0010'] = "rows"
#metaKeys['0028|0011'] = "columns"
#metaKeys['0028|0008'] = "number of frames"
#metaKeys['0054|0081'] = "number of slices"
#metaKeys['0020|1002'] = "images in acquisition"
#metaKeys['0018|0050'] = "slice thickness"
#metaKeys['0028|0030'] = "voxel spacing"
#metaKeys['0008|103e'] = "series description"
#metaKeys['0018|1030'] = "protocol name"
#metaKeys['0008|0070'] = "manufacturer"
#metaKeys['0008|1090'] = "manufacturer's model name"
#metaKeys['0018|1000'] = "device serial number"

def readDICOMMetadataSITK(inputDirectory):
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(inputDirectory)
    reader.SetFileNames(dicom_names)
    reader.MetaDataDictionaryArrayUpdateOn()
    reader.LoadPrivateTagsOn()
    reader.Execute()

    # Get the sorted file names, opens all files in the directory and reads the meta-information
    # without reading the bulk pixel data
    series_ID = reader.GetMetaData(0,'0020|000e')
    sorted_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(inputDirectory, series_ID)
    #print(len(sorted_file_names))

    # Read the bulk pixel data
    #img = sitk.ReadImage(sorted_file_names)
    thisDict = OrderedDict()
    patID = None
    studyDate = None
    studyID = None
    try: 
        patID = reader.GetMetaData(0,metaIDs["patientID"])
    except:
        pass
    try:
        studyDate = reader.GetMetaData(0, metaIDs["studyDate"])
    except:
        pass
    try:
        studyID = reader.GetMetaData(0, metaIDs["studyID"])
    except:
        pass
    
    #thisDict["number of slices"] = len(sorted_file_names)
    for k,v in metaKeys.items():
        try:
            thisDict[v]=reader.GetMetaData(0,k)
        except:
            pass
    return patID, studyDate, studyID, thisDict


def main():
    '''This script receives the path to a folder and recursively finds all dicom files, reporting metadata found'''

    parser = argparse.ArgumentParser()
    parser.add_argument("-idir", help="path to root", default="", required=True)
    #parser.add_argument("-ifilegraph", help="path to graph to convert to skeleton", default=None, required=False)
    parser.add_argument("-ofile", help="name to output json", default=None, required=False)
    parser.add_argument("-special", help="which folder are we segmenting", default=None, required=False)
    parser.add_argument("-contain", help="seek for something in path", default=None, required=False)


    args = parser.parse_args()
    idir = args.idir
    ojson = args.ofile
    special = args.special
    contain = args.contain
    #outputfgraph = args.ofilegraph
    #inputfgraph = args.ifilegraph

    # if inputfgraph is None:
    #     binSketoSke(inputf, outputf, outputfgraph)
    # else:
    #     graphtoSke(inputfgraph, outputf)
    #cases = glob.glob(idir+os.path.sep+"*"+os.path.sep+"PRE"+os.path.sep+"*"+extensionread)
    unknowns = 1
    dirDict = OrderedDict()
    print("Kevin" in special)
    if special is not None:
        print("lvl1")
        if "AP" in special:
            print("ap")
            for root, subd, _ in os.walk(idir): #recursively walk over the directory
                if len(subd)==0 and "PRE" in root:
                    #os.path.join(root, name)
                    galgoID = root.split(idir)[-1].split(os.path.sep)[1]
                    _, _, _, thisDict = readDICOMMetadataSITK(root)
                    dirDict[galgoID] = thisDict #overwrites if multiple series from same patient. fix later
                    print(galgoID)
        elif "Kevin" in special:
            print("kevin")
            for root, subd, files in os.walk(idir): #recursively walk over the directory
                if len(subd)==0 and len(files)!=0:# and "PRE" in root:
                    print(root)
                    #os.path.join(root, name)
                    kevID = root.split(idir)[-1].split(os.path.sep)[1]
                    patID, studyDate, studyID, thisDict = readDICOMMetadataSITK(root)
                    if kevID not in dirDict.keys():
                        dirDict[kevID] = OrderedDict()
                    if patID is not None:
                        if patID not in dirDict.keys():
                            dirDict[kevID][patID] = OrderedDict()
                        dirDict[kevID][patID][studyDate] = thisDict
                    else:
                        if studyID is not None:
                            dirDict[kevID][studyID] = thisDict
                        else:
                            print("something did not work with", root)
                    #dirDict[galgoID] = thisDict #overwrites if multiple series from same patient. fix later
                    #print(galgoID)
    
    
    #    pass
        #print root, d_names, f_names
    if ojson is not None:
        oFile = os.path.join(idir, ojson)
        print(oFile)
        json_dict = OrderedDict()
        #json_dict["name"] = json_name
        json_dict["description"] = "All data in folder "+idir
        timestamp = datetime.today()
        json_dict["timestamp"] = str(timestamp)
        #json_dict["task"] = json_task
        #json_dict["author"] = json_author
        json_dict["studies"] = dirDict
        json_dict["id"] = hashlib.md5(json.dumps(json_dict).encode("utf-8")).hexdigest()[:12]
        with open(oFile, 'w') as f:
            json.dump(json_dict, f, sort_keys=True, indent=4)

    print(dirDict)
    # patID, studyID, thisDict = readDICOMMetadataSITK(idir)

    # if patID is not None:#in dirDict.keys():
    #     if patID not in dirDict.keys():
    #         dirDict[patID] = OrderedDict()
    #     dirDict[patID][studyID] = thisDict
    # else:
    #     if studyID is not None:
    #         dirDict[studyID] = thisDict
    #     else:
    #         dirDict["unknown"+unknowns] = thisDict
    #         unknowns=unknowns+1
    
    # print (dirDict)


    # if ojson is not None:
    #     oFile = os.path.join(idir, ojson)
    #     print(oFile)
    #     json_dict = OrderedDict()
    #     #json_dict["name"] = json_name
    #     json_dict["description"] = "All data in folder "+idir
    #     timestamp = datetime.today()
    #     json_dict["timestamp"] = str(timestamp)
    #     #json_dict["task"] = json_task
    #     #json_dict["author"] = json_author
    #     json_dict["studies"] = dirDict
    #     json_dict["id"] = hashlib.md5(json.dumps(json_dict).encode("utf-8")).hexdigest()[:12]
    #     with open(oFile, 'w') as f:
    #         json.dump(json_dict, f, sort_keys=True, indent=4)



if __name__ == "__main__":
    main()