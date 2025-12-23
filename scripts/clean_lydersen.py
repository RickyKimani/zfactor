import json
import os

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    input_path = os.path.join(project_root, 'data', 'lydersen_raw.json')
    output_path = os.path.join(project_root, 'data', 'lydersen_extracted_vals_raw.json')

    print(f"Reading from {input_path}")
    
    try:
        with open(input_path, 'r') as f:
            raw_data = json.load(f)
    except FileNotFoundError:
        print(f"Error: Could not find {input_path}")
        return

    extracted_data = {}

    if 'datasetColl' not in raw_data:
        print("Error: 'datasetColl' not found in input JSON")
        return

    for dataset in raw_data['datasetColl']:
        name = dataset.get('name', 'unknown')
        points = []
        
        if 'data' in dataset:
            for point in dataset['data']:
                if 'value' in point:
                    points.append(point['value'])
        
        # Sort points by Pr (index 0)
        points.sort(key=lambda x: x[0])
        
        extracted_data[name] = points
        print(f"Extracted {len(points)} points for {name}")

    print(f"Writing to {output_path}")
    with open(output_path, 'w') as f:
        json.dump(extracted_data, f, indent=4)

    print("Done.")

if __name__ == "__main__":
    main()
