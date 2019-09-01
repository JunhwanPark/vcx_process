import sys
import argparse
import os
import json
import xml.etree.ElementTree as ET
import math
import glob

def create_parser():
    parser = argparse.ArgumentParser(description='Calculate VCX Score')
    parser.add_argument('folder', type=str, help='Folder where TE42 xmls are')
    parser.add_argument('--config', type=str, help='Config file', default='vcx1.5_config.json')
    return parser

def read_from_xml(root, xml_entry):
    child = root.find(".//" + xml_entry)
    return child.text

def normalizeToNyquist(mtf, nyquist):
    if mtf > nyquist:
        mtf = nyquist
    elif mtf < 0:
        mtf = nyquist
    return mtf

def get_value_from_xml(xml_entry, xmlFile, valueType):
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

def find_subScoreFolder(folder, testcase):
    folders = [x[0] for x in os.walk(folder)]
    for folder in folders:
        if testcase in folder:
            return folder
    sys.exit("Connot find IQAnalyzer xml for testcase " + testcase)

def calculateScore(value, formula, LGC, HGC, weight):
    if formula == "logarithmic":
        if value < LGC:
            result = 0
        elif value > HGC:
            result = 1
        else:
            result = math.log(value-LGC+1,HGC-LGC+1)
    elif formula == "flat_roof":
        if value > LGC:
            result = 0
        elif value < HGC:
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
    else:
        sys.exit("Formula not supported.")
    result *= weight
    return result

def calculateFinalSubScore(metrics):
    scoreSum = 0
    weightSum = 0
    for name, metric in metrics.items():
        scoreSum += metric[0]
        weightSum += metric[1]
    return (scoreSum/weightSum)*100

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args(sys.argv[1:])
	
    # Find config file
    if not os.path.isfile(args.config):
        sys.exit("Config file does not exist.")
    
    
    with open(args.config, 'r') as f:
        config = json.load(f)

    for subscore in config["IQ_SubScores"]:
        print(subscore["name"])
        
        scores = list()
        
        # Find image xml file based on SubScore name
        subScoreFolder = find_subScoreFolder(args.folder,subscore["name"])
        for xmlFile in glob.glob(subScoreFolder + "\*.xml"):
            print(xmlFile)
            
            results = dict()
            
            for metric in subscore["Metrics"]:
                
                
                # Get calue from xml file
                value = get_value_from_xml(metric["xml_entry"], xmlFile, metric["valueType"])
                
                #value = 21
                
                # Calculate score
                metricScore = calculateScore(value, metric["formula"], metric["LGC"], metric["HGC"], metric["weight"])
                #print(metric["name"], str(value), str(metricScore))
                
                results[metric["name"]] = (metricScore, metric["weight"])
                
            
            scores.append(calculateFinalSubScore(results))
            
        print(subscore["name"], max(scores))
    
    # Perfomance scores
    for metric in config["Performance"]["Metrics"]:
        print(metric["name"])