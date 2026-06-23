import csv
import os
import sys

def csv_to_schema(input_csv):
    output_config = 'config.yml'
    output_schema = 'config.schema.yml'

    try:
        with open(input_csv, mode='r', encoding='utf-8') as csvfile, \
             open(output_schema, mode='w', encoding='utf-8') as schema_file, \
             open(output_config, mode='w', encoding='utf-8') as config_file:

            reader = csv.DictReader(csvfile)
            schema_file.write('$schema: "http://json-schema.org/draft-04/schema#"\n')
            schema_file.write('description: Configuration schema\n')
            schema_file.write('properties:\n')
            names_list = []

            for row in reader:
                name_val = row.get('Name', '').strip()
                type_val = row.get('Type', '').lower()
                default_val = row.get('Default', '').strip('"').strip("'")
                description_val = row.get('Description', '').strip('"').strip("'").replace('"', "'")
                min_val = row.get('Minimum', '').strip()
                max_val = row.get('Maximum', '').strip()
                enum_val = row.get('Enum', '').strip()
                pattern_val = row.get('Pattern', '').strip('"').strip()

                if not name_val:
                    print("Warning: Name field is missing in the CSV.")
                    continue
                if type_val.lower() not in ["integer", "boolean", "number", "string"]:
                    print('Error: Type value is not ["integer", "boolean", "number", "string"]')
                    sys.exit(1)

                names_list.append(name_val)

                # Write to schema file
                schema_file.write(f"  {name_val}:\n")
                schema_file.write(f"    type: {type_val}\n")
                schema_file.write(f"    description: \"{description_val}\"\n")

                # Handle integers and numbers
                if type_val.lower() in ["integer", "number"]:
                    config_file.write(f"{name_val}: {default_val}\n")
                    schema_file.write(f"    default: {default_val}\n")
                    if min_val:
                        schema_file.write(f"    minimum: {min_val}\n")
                    if max_val:
                        schema_file.write(f"    maximum: {max_val}\n")
                        
                # Handle booleans
                elif type_val.lower() == "boolean":
                    config_file.write(f"{name_val}: {default_val.lower().capitalize()}\n")
                    schema_file.write(f"    default: {default_val.lower().capitalize()}\n")

                # Handle strings        
                elif type_val.lower() == "string":
                    config_file.write(f"{name_val}: \"{default_val}\"\n")
                    schema_file.write(f'    default: \"{default_val}\"\n')
                    if enum_val:
                        schema_file.write(f"    enum: {enum_val}\n")
                    if pattern_val:
                        schema_file.write(f"    pattern: \"{pattern_val}\"\n")
                
                schema_file.write("\n")
            
            schema_file.write("required:\n")
            for name_val in names_list:
                schema_file.write(f"  - {name_val}\n")

        print(f"Successfully created {output_config} and {output_schema} in the current directory.")

    except FileNotFoundError:
        print(f"Error: The file {input_csv} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print("Usage: python script.py <input_csv_path>\nConvert a config.schema.csv file into a config.yml and config.schema.yml in the same directory as the python script is running.")
        sys.exit(1)

    input_csv_path = sys.argv[1]
    csv_to_schema(input_csv_path)
