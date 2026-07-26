"""Convert a JSON Schema-style config.schema.yml into a YAML default config.yml.

This script reads the top-level `properties` section from a schema file and
writes a config file containing each property name with its `default` value.
"""

import argparse
import sys

try:
    import yaml
except ImportError as exc:
    raise SystemExit(
        "PyYAML is required to run this script. Install it with `pip install pyyaml`."
    ) from exc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert config.schema.yml into config.yml format."
    )
    parser.add_argument("schema_path", help="Path to the input config.schema.yml file.")
    parser.add_argument(
        "--output",
        "-o",
        dest="output_path",
        default="config.yml",
        help="Path to write the generated config.yml file (default: config.yml).",
    )
    return parser.parse_args()


def load_schema(path):
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError("Schema file did not contain a valid YAML mapping.")
    return data


class QuotedString(str):
    pass


def quoted_string_representer(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="'")


def bool_representer(dumper, data):
    value = "True" if data else "False"
    return dumper.represent_scalar("tag:yaml.org,2002:bool", value)


yaml.SafeDumper.add_representer(QuotedString, quoted_string_representer)

yaml.SafeDumper.add_representer(bool, bool_representer)


def schema_to_config(schema):
    properties = schema.get("properties")
    if not isinstance(properties, dict):
        raise ValueError("Schema file must include a top-level 'properties' mapping.")

    config = {}
    for key, value in properties.items():
        if not isinstance(value, dict):
            continue

        default = value.get("default")
        if value.get("type") == "string":
            config[key] = QuotedString(default if default is not None else "")
        else:
            config[key] = default
    return config


def write_config(config, output_path):
    with open(output_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(
            config,
            handle,
            sort_keys=False,
            default_flow_style=False,
            allow_unicode=True,
            width=4096,
        )


def main():
    args = parse_args()
    schema = load_schema(args.schema_path)
    config = schema_to_config(schema)
    write_config(config, args.output_path)
    print(f"Generated {args.output_path} from {args.schema_path}.")


if __name__ == "__main__":
    main()
