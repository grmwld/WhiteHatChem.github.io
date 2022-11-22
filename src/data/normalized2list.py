# transform the normalized file into a list of file for each molecule
import json
import os
def generate_json_from_normalized(input_path='normalized.json', output_path='json/'):
    with open(input_path, 'r') as fi:
        data = json.load(fi)
    n = len(data['name'])
    output = {k: {} for k, v in data['name'].items()}
    
    for k1, v1 in data.items():
        for k2, v2 in v1.items():
            output[k2][k1] = v2

    for k, v in output.items():
        f_name = f"{v['name']}.json"
        with open(os.path.join(output_path, f_name), 'w+') as fo:
            json.dump(v, fo, indent=2)



if __name__ == "__main__":
    generate_json_from_normalized()
