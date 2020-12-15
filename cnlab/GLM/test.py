import os
if __name__ == "__main__":

    if os.environ.get('SINGULARITY_COMMAND'):
        inSingularity = True
    print(os.environ)