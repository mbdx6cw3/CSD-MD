import yaml, argparse, os

def csdMD():
    '''

    :returns md_params: contains all molecular dynamics parameters required to
                        run a simulation.
    '''
    # parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--md_params", default="input.yaml")
    args = parser.parse_args()
    input_file = args.md_params
    isExist = os.path.exists(input_file)
    if not isExist:
        print("Error: input file does not exist.")
    else:
        with open(input_file, 'r') as input:
            md_params = yaml.full_load(input)

    return md_params

