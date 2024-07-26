import json
import argparse
import os
import csv
from angiographies.utils.fileoperations import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ifolder", type=str, required=True, help="Folder to read")
    #parser.add_argument("-date", type=str, required=False, default="None",help="csv date")
    #parser.add_argument("-ignore", type=str, required=False, help="use this if you want to run on 2d images", default="")
    parser.add_argument("--overunder", help="compute over under estimation", action="store_true", required=False, default=False)
    parser.add_argument("--size", help="compute over nidus size", action="store_true", required=False, default=False)
    #parser.add_argument("--cuda", help="use this if you want to try running on gpu", action="store_true", required=False, default=False)
    
    args = parser.parse_args()
    ifolder = args.ifolder
    #append = args.a
    segm1=["dl"]#, "manual"]#,"dldice"]
    segm2=["dl_cldice"]
    translate = {"morphological":"Morph", "skeletonedt": "Thin", "vmtk": "VMTK", "boundingbox" : "BB", "hull":"CH", "spheres" : "Sph"}

    header = ["Case", "Segmentation","Metric", "Value","Method", "Tries"] #header manual
    if args.size:
        header = ["Case", "Segmentation","Size Cat GT", "Size Cat Pred","Size GT", "Size Pred","Overestimation", "Underestimation","Method", "Tries"] #header manual
    alldices = []
    for segm in segm2:
        allmethods = getsubdirs(ifolder, contains = segm+"_", join=False) if isdir(ifolder) else None
        for evaluation in allmethods:
            method = evaluation.split("_")[2] #1 if doing five cases, 2 if seventeen
            method = translate[method] if method in translate else method
            submethod = evaluation.split("_")[3] #2 if doing five cases, 3 if seventeen
            submethod = translate[submethod] if submethod in translate else submethod
            with open(os.path.join(ifolder, evaluation,"summary.json"),"r") as f:
            
                # returns JSON object as
                # a dictionary
                data = json.load(f)
                
                # Iterating through the json
                # list
                for i in data['results']['all']:
                    intentos = 0
                    case = (i['reference'].split(os.path.sep)[-1]).split(".nii.gz")[0]
                    if (case in ["3009", "3034"] and segm == "dl" and method == "Thin") or (case in ["3020"] and segm == "manual" and method == "Thin") or \
                        (case in ["3001", "3033"] and segm == "dl_cldice" and method == "Thin"):
                        intentos = 1
                    if (case in ["3032"] and segm == "manual" and method == "Thin"):
                        intentos = 3
                    if (case in ["3002"] and segm == "dl_cldice" and method == "Thin"):
                        intentos = 5
                    if (case in ["3003"] and segm == "dl_cldice" and method == "Thin"):
                        intentos = 2
                    if args.size:
                        sizegt=int(i["1"]["Total Positives Reference"])*(0.0227**3)
                        sizepred=int(i["1"]["Total Positives Test"])*(0.0227**3)
                        if (sizegt)<14.1371:
                            gradeGT = 1
                        elif (sizegt)>113.0973:
                            gradeGT = 3
                        else:
                            gradeGT = 2
                        if (sizepred)<14.1371:
                            gradePred = 1
                        elif (sizepred)>113.0973:
                            gradePred = 3
                        else:
                            gradePred = 2
                        alldices.append([case, segm, gradeGT, gradePred,sizegt, sizepred,float(i["1"]["False Positives"])/float(i["1"]["Total Positives Reference"]), \
                                        float(i["1"]["False Negatives"])/float(i["1"]["Total Positives Reference"]), method + "-" + submethod, intentos])
                    else:
                        for metric in ["Dice", "Recall", "Precision"]:
                            value = float(i["1"][metric])
                            alldices.append([case, segm, metric, value,method + "-" + submethod, intentos])
                        if args.overunder:
                            alldices.append([case, segm, "Overestimation", float(i["1"]["False Positives"])/float(i["1"]["Total Positives Reference"]), gradeGT, gradePred,method + "-" + submethod, intentos])
                            alldices.append([case, segm, "Underestimation", float(i["1"]["False Negatives"])/float(i["1"]["Total Positives Reference"]), gradeGT, gradePred,method + "-" + submethod, intentos])
                            if (int(i["1"]["Total Positives Reference"])*(0.0227**3))<14.1371:
                                grade = 1
                            elif (int(i["1"]["Total Positives Reference"])*(0.0227**3))>113.0973:
                                grade = 3
                            else:
                                grade = 2
                            alldices.append([case, segm, "Size GT", grade, method + "-" + submethod, intentos])
                            if (int(i["1"]["Total Positives Test"])*(0.0227**3))<14.1371:
                                grade = 1
                            elif (int(i["1"]["Total Positives Test"])*(0.0227**3))>113.0973:
                                grade = 3
                            else:
                                grade = 2
                            alldices.append([case, segm, "Size Pred", grade, method + "-" + submethod, intentos])

                        
    print(header)
    print(alldices)

    with open(os.path.join(ifolder,"svtdices-fixed-morph-size.csv"), 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
        # write multiple rows
        writer.writerows(alldices)



if __name__ == "__main__":
    main()