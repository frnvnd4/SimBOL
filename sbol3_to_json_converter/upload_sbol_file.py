from google.colab import files

def upload_sbol_file():
    print("Please select your SBOL 3 file to upload:\n")

    uploaded = files.upload()

    if uploaded:
        sbol_file = list(uploaded.keys())[0]
        name = sbol_file.rsplit(".", 1)[0]
        print(f"------- File uploaded successfully: {sbol_file} -------")
        return sbol_file, name
    else:
        raise Exception("------- No file was uploaded. Please try again. -------")
