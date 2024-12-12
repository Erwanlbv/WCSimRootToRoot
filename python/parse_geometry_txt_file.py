#!/usr/bin/env python3

import sys
import yaml

def parse_text_file(file_path):
        
    data_dict = {
        "Detector": {},
        "Type 1 PMT": {},
        "Type 2 PMT": {},
        "OD PMT": {},
        "Centre offset": [],
        "Entries": []
    }

    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        for line in lines:

            # The first part deals with the heading of the .txt file
            if "Detector radius & height" in line:
                radius, height = line.split()[-2:]
                data_dict["Detector"]["Radius"] = float(radius)
                data_dict["Detector"]["Height"] = float(height)

            
            elif "Type 1 ID PMT number & size" in line:
                number, size = line.split()[-2:]
                data_dict["Type 1 PMT"]["Number"] = int(number)
                data_dict["Type 1 PMT"]["Size"] = float(size)

            
            elif "Type 2 ID PMT number & size" in line:
                number, size = line.split()[-2:]
                data_dict["Type 2 PMT"]["Number"] = int(number)
                data_dict["Type 2 PMT"]["Size"] = float(size)
            
            elif "OD PMT number & size" in line:
                _, number, size = line.rsplit("  ", 2)
                data_dict["OD PMT"]["Number"] = int(number)
                data_dict["OD PMT"]["Size"] = float(size)
            
            elif "Centre offset" in line:
                _, values = line.split("  ", 1)
                offsets = [float(val) for val in values.split()]
                data_dict["Centre offset"] = offsets

            # The second part deals with the body (the PMTs informations) of the .txt file
            else:
                parts = line.split()
                if len(parts) == 10: # Maybe HK-only ?
                    entry = {
                        "id": int(parts[0]),
                        "sub_id": int(parts[1]),
                        "mPMT_number": int(parts[2]),
                        "position": [float(parts[3]), float(parts[4]), float(parts[5])]
                    }
                    data_dict["Entries"].append(entry)
    return data_dict


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_file.py path_to_txt_file.txt path_to_save_output.yaml")
        sys.exit(1)
    
    file_path, output_path = sys.argv[1], sys.argv[2]
    print(f"\n\nParsing data from {file_path}")
    data_dict = parse_text_file(file_path)
    
    # Serialize data_dict to YAML and write it to a file
    with open(output_path, 'w') as f:
        yaml.dump(data_dict, f, default_flow_style=False)
    
    print(f"Data saved at {output_path}\n\n")
