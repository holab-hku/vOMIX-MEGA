import csv
import os
import sys

def csv_to_schema(input_csv):
    output_config = 'config.yml'
    output_schema = 'config.schema.yml'

    try:
        # Using utf-8-sig automatically maps BOM sequences cleanly
        with open(input_csv, mode='r', encoding='utf-8-sig') as csvfile, \
             open(output_schema, mode='w', encoding='utf-8') as schema_file, \
             open(output_config, mode='w', encoding='utf-8') as config_file:

            reader = csv.DictReader(csvfile)
            
            # Clean up fieldnames just in case there's unexpected whitespace
            if reader.fieldnames:
                reader.fieldnames = [field.strip() for field in reader.fieldnames]

            schema_file.write('$schema: "http://json-schema.org/draft-04/schema#"\n')
            schema_file.write('description: Configuration schema\n')
            schema_file.write('properties:\n')
            names_list = []

            for row in reader:
                name_val = row.get('Name', '').strip()
                type_val = row.get('Type', '').lower().strip()
                
                # Keep raw strings intact without stripping inner target characters
                default_val = row.get('Default', '')
                description_val = row.get('Description', '').replace('"', "'")
                min_val = row.get('Minimum', '').strip()
                max_val = row.get('Maximum', '').strip()
                enum_val = row.get('Enum', '').strip()
                pattern_val = row.get('Pattern', '').strip()

                if not name_val:
                    continue
                if type_val not in ["integer", "boolean", "number", "string"]:
                    print(f'Error: Type value "{type_val}" is not integer, boolean, number, or string')
                    sys.exit(1)

                names_list.append(name_val)

                # Write to schema file
                schema_file.write(f"  {name_val}:\n")
                schema_file.write(f"    type: {type_val}\n")
                schema_file.write(f"    description: \"{description_val}\"\n")

                # Handle integers and numbers
                if type_val in ["integer", "number"]:
                    config_file.write(f"{name_val}: {default_val.strip()}\n")
                    schema_file.write(f"    default: {default_val.strip()}\n")
                    if min_val:
                        schema_file.write(f"    minimum: {min_val}\n")
                    if max_val:
                        schema_file.write(f"    maximum: {max_val}\n")
                        
                # Handle booleans
                elif type_val == "boolean":
                    clean_bool = default_val.strip().lower().capitalize()
                    config_file.write(f"{name_val}: {clean_bool}\n")
                    schema_file.write(f"    default: {clean_bool}\n")

                # Handle strings safely keeping values as raw characters       
                elif type_val == "string":
                    config_file.write(f"{name_val}: {default_val}\n")
                    schema_file.write(f'    default: {default_val}\n')
                    if enum_val:
                        schema_file.write(f"    enum: {enum_val}\n")
                    if pattern_val:
                        schema_file.write(f"    pattern: {pattern_val}\n")
                
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
    if len(sys.argv) != 2 or sys.argv[1] in ["-h", "--help"]:
        print("Usage: python script.py <input_csv_path>\nConvert a config.schema.csv file into a config.yml and config.schema.yml.")
        sys.exit(1)

    input_csv_path = sys.argv[1]
    csv_to_schema(input_csv_path)
