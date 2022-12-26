import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle

# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

def desc_calc2():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output_cox2.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule_cox2.smi')


# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href
def filedownload2(df2):
    csv = df2.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction_cox2.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('COX_1_FILE.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output for cox1**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)
def build_model2(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('COX_2_FILE.pkl', 'rb'))
    # Apply model to make predictions
    prediction2 = load_model.predict(input_data)
    st.header('**Prediction output for cox2**')
    prediction_output2 = pd.Series(prediction2, name='pIC50')
    molecule_name2 = pd.Series(load_data2[1], name='molecule_name')
    df2 = pd.concat([molecule_name2, prediction_output2], axis=1)
    st.write(df2)
    st.markdown(filedownload2(df2), unsafe_allow_html=True)


# Logo image
image = Image.open('logo.png')

st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Bioactivity Prediction App (Acetylcholinesterase)
This app allows you to predict the bioactivity towards inhibting the `Acetylcholinesterase` enzyme. `Acetylcholinesterase` is a drug target for Alzheimer's disease.
**Credits**
- App built in `Python` + `Streamlit` by [Chanin Nantasenamat](https://medium.com/@chanin.nantasenamat) (aka [Data Professor](http://youtube.com/dataprofessor))
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
---
""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
""")

if st.sidebar.button('Predict_cox1'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors **')
    desc = pd.read_csv('descriptors_output.csv')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
elif st.sidebar.button('Predict_cox2'):
    load_data2 = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data2.to_csv('molecule_cox2.smi', sep = '\t', header = False, index = False)

    st.header('**Original input data**')
    st.write(load_data2)

    with st.spinner("Calculating descriptors for cox2..."):
        desc_calc2()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors for cox2**')
    desc2 = pd.read_csv('descriptors_output_cox2.csv')
    st.write(desc2)
    st.write(desc2.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_cox2_list.csv').columns)
    desc_subset2 = desc2[Xlist]
    st.write(desc_subset2)
    st.write(desc_subset2.shape)

    # Apply trained model to make prediction on query compounds
    build_model2(desc_subset2)

else:
    st.info('Upload input data in the sidebar to start!')