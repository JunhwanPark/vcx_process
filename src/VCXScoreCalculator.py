import sys
import argparse
import os
import json
import xml.etree.ElementTree as ET
import math
import glob
import csv
import numpy as np
import matplotlib.pyplot as plt

def create_parser():
    parser = argparse.ArgumentParser(description='Calculate VCX Score')
    parser.add_argument('folder', type=str, help='Folder where TE42 xmls are')
    parser.add_argument('--config', type=str, help='Config file', default='vcx1.5_config.json')
    parser.add_argument('--plotfunctions', action='store_true', help='Plot functions of all metrics and save them as png files.')
    return parser

def read_from_xml(root, xml_entry):
    child = root.find(".//" + xml_entry)
    if child is None:
        raise Exception("XML entry not found.")
    return child.text

def normalizeToNyquist(mtf, nyquist):
    if mtf > nyquist:
        mtf = nyquist
    elif mtf < 0:
        mtf = nyquist
    return mtf

def get_value_from_xml(xml_entry, xmlFile, valueType):
    
    try:
        root = ET.parse(xmlFile)
        
        # Simple floating point value
        if valueType == "float":
            return float(read_from_xml(root, xml_entry))
        # Maximum of table
        elif valueType == "max":
            data = read_from_xml(root, xml_entry).split(";")
            data = list(filter(None, data))
            data = [float(a) for a in data]
            return max(data)
        # Mean of multiple entries from single xml
        elif valueType == "multiXMLEntryMean":
            data = list()
            for entry in xml_entry:
                data.append(float(read_from_xml(root, entry)))
            return sum(data)/len(data)
        # Calculate means of entries except last. Normalize result with last entry.
        elif valueType == "overshoot":
            data = list()
            for entry in xml_entry[:-1]:
                data.append(float(read_from_xml(root, entry)))
            result = sum(data)/len(data)
            result = (100*result)/float(read_from_xml(root, xml_entry[-1]))
            return result
        # Effective pixel count
        # Calculated from siemens star MTF10 and pixel count
        elif valueType == "EPC":
            nyquist = float(read_from_xml(root, xml_entry[-2]))
            pixelCount = float(read_from_xml(root, xml_entry[-1]))
            center = list()
            for i in range(8):
                value = float(read_from_xml(root, xml_entry[i]))
                # Normalize orthogonal sectors to nyquist and diagonal to 1.4 * nyquist
                if i % 2 == 0:
                    value = normalizeToNyquist(value, nyquist)
                else:
                    value = normalizeToNyquist(value, nyquist*1.4)
                center.append(value)
            corner = list()
            #Corner 1
            for i in range(3):
                value = float(read_from_xml(root, xml_entry[i+8]))
                # Normalize orthogonal sectors to nyquist and diagonal to 1.4 * nyquist
                if i % 2 == 0:
                    value = normalizeToNyquist(value, nyquist)
                else:
                    value = normalizeToNyquist(value, nyquist*1.4)
                corner.append(value)
            #Corner 2
            for i in range(3):
                value = float(read_from_xml(root, xml_entry[i+11]))
                # Normalize orthogonal sectors to nyquist and diagonal to 1.4 * nyquist
                if i % 2 == 0:
                    value = normalizeToNyquist(value, nyquist)
                else:
                    value = normalizeToNyquist(value, nyquist*1.4)
                corner.append(value)
            #Corner 3
            for i in range(3):
                value = float(read_from_xml(root, xml_entry[i+14]))
                # Normalize orthogonal sectors to nyquist and diagonal to 1.4 * nyquist
                if i % 2 == 0:
                    value = normalizeToNyquist(value, nyquist)
                else:
                    value = normalizeToNyquist(value, nyquist*1.4)
                corner.append(value)
            #Corner 4
            for i in range(3):
                value = float(read_from_xml(root, xml_entry[i+17]))
                # Normalize orthogonal sectors to nyquist and diagonal to 1.4 * nyquist
                if i % 2 == 0:
                    value = normalizeToNyquist(value, nyquist)
                else:
                    value = normalizeToNyquist(value, nyquist*1.4)
                corner.append(value)
            
            center = sum(center)/len(center)
            corner = sum(corner)/len(corner)
            
            return (((center+corner)/(2*nyquist))**2)*pixelCount
        # Artifact metric which is calculated from dead leaves chart
        elif valueType == "artifacts":
            DL_cross = float(read_from_xml(root, xml_entry[0]))/100
            DL_direct_old = float(read_from_xml(root, xml_entry[1]))/100
            artifacts = 100-(DL_cross/(DL_direct_old/100))
            if artifacts < 0:
                return 0
            else:
                return artifacts
        # Floating point value, which is normalized to nyquist before reporting
        elif valueType == "normalizeNyquist":
            value = float(read_from_xml(root, xml_entry[0]))
            nyquist = float(read_from_xml(root, xml_entry[1]))
            return normalizeToNyquist(value, nyquist)
        # Mean value of certain skin tone color patches
        elif valueType == "skinMean":
            patchNos = [4,9,62,63,64,65,66,67,68,74,75,76,77,78,79,80,92]
            patches = read_from_xml(root, xml_entry).split(";")
            
            mean = 0.0
            for number in patchNos:
                mean += float(patches[number])
            mean /= len(patchNos)
            
            return mean
        # Only for debugging purposes
        elif valueType == "debug":
            return 2.2
        else:
            sys.exit("Unknown value type.")
    except:
        print("XML entry not found")

def find_subScoreFolder(folder, testcase):
    folders = [x[0] for x in os.walk(folder)]
    for folder in folders:
        if testcase in folder:
            return folder
    print("Cannot find IQAnalyzer xml for testcase " + testcase)
    return ""

def calculateScoreV20(value, formula, LGC, HGC, weight):
    if formula == "logarithmic":
        if value < LGC:
            result = 0
        elif value > HGC:
            result = 1
        else:
            result = math.log(value-LGC+1,HGC-LGC+1)
    elif formula == "flat_roof":
        if abs(value) > LGC:
            result = math.log(abs(value-LGC),LGC)*(-2)
            return result
        elif abs(value) < HGC:
            result = 1
        else:
            result = (LGC-abs(value))/(LGC-HGC)
    elif formula == "logarithmic_neg.linear":
        if value < 2*LGC-HGC:
            result = -1
        elif value < LGC:
            result = (LGC-value)/(LGC-HGC)
        elif value > HGC:
            result = 1
        else:
            result = math.log(value-LGC+1,HGC-LGC+1)
    elif formula == "roof_negative_ll":
        if value < 2*LGC-HGC or value > 3*HGC - 2*LGC:
            result = -1
        else:
            result = 1-abs(value-HGC)/(HGC-LGC)
    elif formula == "linear":
        if value < LGC:
            result = math.log(abs(value-LGC),LGC)*(-1)
            return result
        elif value > HGC:
            result = 1
        else:
            result = (LGC-value)/(LGC-HGC)
    elif formula == "roof_hl":
        if value < 2*HGC-LGC or value > LGC:
            result = math.log(abs(value-LGC),LGC)*(-1)
            return result
        else:
            result = 1-abs(value-HGC)/(LGC-HGC)
    elif formula == "logarithmic_roof":
        if abs(value) > LGC:
            result = math.log(abs(value-LGC),LGC)*(-2)
        elif abs(value) < HGC:
            result = 1
        else:
            result = 1 - math.log(abs(value)-HGC+1,LGC-HGC+1)
    elif formula == "roof_II":
        if value < LGC or value > 2*HGC-LGC:
            result = abs(value-LGC)*(-2)
            return result
        else:
            result = 1-abs(value-HGC)/(HGC-LGC)
    else:
        sys.exit("Formula not supported.")
    result *= weight
    return result

def calculateScoreV15(value, formula, LGC, HGC, weight):
    if formula == "logarithmic":
        if value < LGC:
            result = 0
        elif value > HGC:
            result = 1
        else:
            result = math.log(value-LGC+1,HGC-LGC+1)
    elif formula == "flat_roof":
        if abs(value) > LGC:
            result = 0
        elif abs(value) < HGC:
            result = 1
        else:
            result = (LGC-abs(value))/(LGC-HGC)
    elif formula == "logarithmic_neg.linear":
        if value < 2*LGC-HGC:
            result = -1
        elif value < LGC:
            result = (LGC-value)/(LGC-HGC)
        elif value > HGC:
            result = 1
        else:
            result = math.log(value-LGC+1,HGC-LGC+1)
    elif formula == "roof_negative_ll":
        if value < 2*LGC-HGC or value > 3*HGC - 2*LGC:
            result = -1
        else:
            result = 1-abs(value-HGC)/(HGC-LGC)
    elif formula == "linear":
        if value < LGC:
            result = 0
        elif value > HGC:
            result = 1
        else:
            result = (LGC-value)/(LGC-HGC)
    elif formula == "roof_hl":
        if value < 2*HGC-LGC or value > LGC:
            result = 0
        else:
            result = 1-abs(value-HGC)/(LGC-HGC)
    elif formula == "logarithmic_roof":
        if abs(value) > LGC:
            result = 0
        elif abs(value) < HGC:
            result = 1
        else:
            result = 1 - math.log(abs(value)-HGC+1,LGC-HGC+1)
    elif formula == "roof_II":
        if value < LGC or value > 2*HGC-LGC:
            result = 0
        else:
            result = 1-abs(value-HGC)/(HGC-LGC)
    else:
        sys.exit("Formula not supported.")
    result *= weight
    return result

def calculateScore(value, xp, yp, weight):
    score = np.interp(value, xp, yp)
    return score*weight

def calculateFinalSubScore(metrics):
    scoreSum = 0
    weightSum = 0
    for name, metric in metrics.items():
        scoreSum += metric[0]
        weightSum += metric[1]
    return (scoreSum/weightSum)*100

def calculateMeanFileSize(folder):
    files = glob.glob(folder + "\\*.jpg")
    fileSizes = list()
    
    Keywords = ["check", "calc", "test", "analysis"]
    
    # Discard files that include certain keywords. These are IQAnalyzer analysis files.
    files = [s for s in files if not any(k in s for k in Keywords)]
    
    for file in files:
        fileSizes.append(os.path.getsize(file)/1024)
    return sum(fileSizes)/len(fileSizes)

def calculatePerformanceMetric(metric, folders):
    if metric["valueType"] == "compressionLoss":
        framerateFolder = args.folder + "\\" + folders["framerate"]
        AF800lxFolder = args.folder + "\\" + folders["AF_800lux"]
        nonAFFolder = args.folder + "\\" + folders["nonAF"]
        StartupFolder = args.folder + "\\" + folders["startup"]
        
        framerateImageSize = calculateMeanFileSize(framerateFolder)
        AFImageSize = calculateMeanFileSize(AF800lxFolder)
        nonAFImageSize = calculateMeanFileSize(nonAFFolder)
        startupImageSize = calculateMeanFileSize(StartupFolder)
        
        return (framerateImageSize/(sum([AFImageSize, nonAFImageSize, startupImageSize])/3)-1)*100
    
    if metric["valueType"] == "afFailure":
        folder = args.folder + "\\" + folders[metric["fileTag"]]
        
        xmls = glob.glob(folder + "\*.xml")
        
        MTFs = list()
        for xml in xmls:
            meanMTF = 0
            for edge in ["Top", "Bottom", "Left", "Right"]:
                meanMTF += get_value_from_xml("SFR_Edge/EdgesCenter/" + edge + "/Y/Results/Limit/MTF50", xml, "float")
            MTFs.append(meanMTF/4)
            
        # Get nyquist from first xml
        nyquist = get_value_from_xml("SFR_Edge/EdgesCenter/Y/SFR/nyquist_frequency", xml, "float")
        
        # Threshold for blurriness
        threshold = 0.5
        
        unsharpCount = 0
        for MTF in MTFs:
            if MTF < threshold*nyquist:
                unsharpCount += 1
            
        return (unsharpCount/len(MTFs))*100
        
    if metric["valueType"] == "delta":
        # Read result csv
        file = glob.glob(args.folder + "\*" + metric["fileTag"] + "*.csv")[0]
        values = list()
        values2 = list()
        with open(file) as csvFile:
            readCSV = csv.reader(csvFile, delimiter=',')
            for row in readCSV:
                try: 
                    values.append(float(row[metric["column"][0]]))
                    values2.append(float(row[metric["column"][1]]))
                except ValueError:
                    continue
                
        return sum(values)/len(values)-sum(values2)/len(values2)
    
    if metric["valueType"] == "MultipleXMLMean":
        
        value = list()
        for xml in glob.glob(args.folder + "\\" + folders[metric["fileTag"]] + "\\*.xml"):            
            imageMeanValue = list()
            # Read xmls
            for xmlEntry in metric["xml_entry"]:
                imageMeanValue.append(get_value_from_xml(xmlEntry, xml, "float"))
                
            imageMeanValue = sum(imageMeanValue)/len(imageMeanValue)
            
            value.append(imageMeanValue)
            
        value = sum(value)/len(value)
        return(value)
        
    if "MultipleXMLMeanDelta" in metric["valueType"]:
        
        value1 = list()
        for xml in glob.glob(args.folder + "\\" + folders[metric["fileTag"][0]] + "\\*.xml"):            
            imageMeanValue = list()
            # Read xmls
            for xmlEntry in metric["xml_entry"]:
                imageMeanValue.append(get_value_from_xml(xmlEntry, xml, "float"))
                
            imageMeanValue = sum(imageMeanValue)/len(imageMeanValue)
            
            value1.append(imageMeanValue)
            
        value1 = sum(value1)/len(value1)
        
        value2 = list()
        for xml in glob.glob(args.folder + "\\" + folders[metric["fileTag"][1]] + "\\*.xml"):            
            imageMeanValue = list()
            # Read xmls
            for xmlEntry in metric["xml_entry"]:
                imageMeanValue.append(get_value_from_xml(xmlEntry, xml, "float"))
                
            imageMeanValue = sum(imageMeanValue)/len(imageMeanValue)
            
            value2.append(imageMeanValue)
            
        value2 = sum(value2)/len(value2)
        
        if "Percent" in metric["valueType"]:
            return((value1/value2-1)*100)
        if "Absolute" in metric["valueType"]:
            return(value1-value2)
        sys.exit("Value type not supported.")
    
    # Read result csv
    file = glob.glob(args.folder + "\*" + metric["fileTag"] + "*.csv")[0]
    values = list()
    with open(file) as csvFile:
        readCSV = csv.reader(csvFile, delimiter=',')
        for row in readCSV:
            try: 
                values.append(float(row[metric["column"]]))
            except ValueError:
                continue
    
    if metric["valueType"] == "fps":
        return (values[-1] - values[1])/values[0]
    if metric["valueType"] == "meanScaled":
        return (sum(values[1:])/(len(values)-1))/values[0]
    if metric["valueType"] == "mean":
        return sum(values)/len(values)
    sys.exit("Value type not supported.")

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args(sys.argv[1:])
    
    if args.plotfunctions:
        print("Saving plots of metric to score functions.")
	
    # Find config file
    if not os.path.isfile(args.config):
        sys.exit("Config file does not exist.")
    
    
    with open(args.config, 'r') as f:
        config = json.load(f)

    for subscore in config["IQ_SubScores"]:
        print(subscore["name"])
        
        scores = list()
        
        subScoreWeight = subscore["SubScoreWeight"]
        
        # Find image xml file based on SubScore name
        subScoreFolder = find_subScoreFolder(args.folder,subscore["name"])
        
        if(subScoreFolder == ""):
            continue
        
        xmlFiles = glob.glob(os.path.join(subScoreFolder, "*.xml"))
        
        if len(xmlFiles) == 0:
            print("Cannot find IQAnalyzer xml for testcase " + subscore["name"])
            continue
        
        for xmlFile in xmlFiles:
            print(xmlFile)
            
            results = dict()
            
            for group in subscore["Groups"]:
                groupWeight = group["groupWeight"]
                for metric in group["Metrics"]:
                    
                    if args.plotfunctions:
                        plt.clf()
                        plt.plot(metric["xp"], metric["yp"])
                        plt.title(subscore["name"] + " - " + metric["name"])
                        plt.xlabel("Metric")
                        plt.ylabel("Score")
                        plt.savefig(subscore["name"] + " " + metric["name"] + ".png")
                    
                    # Get calue from xml file
                    value = get_value_from_xml(metric["xml_entry"], xmlFile, metric["valueType"])
                    
                    #value = 21
                    
                    finalWeight = subScoreWeight*groupWeight*metric["weight"]*10
                    
                    # Calculate score
                    #metricScore = calculateScoreV20(value, metric["formula"], metric["LGC"], metric["HGC"], finalWeight)
                    metricScore = calculateScore(value, metric["xp"], metric["yp"], finalWeight)
                    if math.isnan(metricScore):
                        metricScore = 0
                    
                    #print(metric["name"], str(value), str(metricScore), str(finalWeight))                
                    results[metric["name"]] = (metricScore, finalWeight)
            
                
            scores.append(calculateFinalSubScore(results))
        
        
        print(subscore["name"], max(scores))
    
    perfMetrics = dict()
    # Perfomance scores
    #print("Performance")
    #for metric in config["Performance"]["Metrics"]:
    #    value = calculatePerformanceMetric(metric, config["Performance"]["Folders"])
    #    metricScore = calculateScoreV20(value, metric["formula"], metric["LGC"], metric["HGC"], metric["weight"])
    #    perfMetrics[metric["name"]] = (metricScore, metric["weight"])
        
    #print("Performance", calculateFinalSubScore(perfMetrics))