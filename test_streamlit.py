# import the necessary packages
import streamlit as st
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
import streamlit.components.v1 as components
import seaborn as sns

# set the web app page name (optional)
st.set_page_config(page_title='SAMPLE WEB APP', page_icon=None, layout="wide")

# set markdown
st.markdown('# First Streamlit UI Web Application by Shilpa.')

data = pd.read_csv ('Housing.csv')
data.columns = data.columns.str.upper() # convert the columns to uppercase

# write a text or comment ( the number of "#" symbols used determines the text size.)
st.write( '### 1. Overview of the Data ')

# view the dataframe in streamlit
st.dataframe(data, use_container_width=True)

st.write( '### 2. Understanding the Data ')

# creating radio button
selected = st.radio( "**What do you want to know about the data ?**", 
["Description", "Data Sample", "Data Head/Tail","Data Shape"])

if selected == 'Description':
    st.dataframe(data.describe(), use_container_width=True)  # shows the basic data descriptive
elif selected == 'Data Sample':
    st.dataframe(data.sample(10), use_container_width=True)  # select random rows

elif selected == 'Data Head/Tail':
    st.dataframe(data.head(), use_container_width=True)  # shows the head of the dataframe 
else:
    st.write('###### The shape of the data is :',data.shape)  # shows the data shape


st.write( '### 3. Data Visualization ')
############# Option-1 : embed and display the html analysis result already done by the auto-eda libraries(pandasprofile, sweet viz) using streamlit Tabs ###############################

# creating tabs 


numerical_features = data.select_dtypes(exclude= 'object').columns 
categorical_features = data.select_dtypes('object').columns 

# creating mutiselect tab in the left sidebar
graph_options = st.sidebar.multiselect( "SELECT GRAPH TYPE:", options=['HISTOGRAM','COUNT-PLOT','BOX-PLOT'], default=['HISTOGRAM','COUNT-PLOT','BOX-PLOT'])

# creating three columns to display on the streamlit web 
col1,col2,col3 = st.columns(3)
# defining columns to disply in the streamlit
col1.write( '##### Histogram Plot ')
col2.write( '##### Box Plot ')
col3.write( '##### Count Plot')

# if histogram selected show the histo graph for reach numerical feature
if  'HISTOGRAM' in graph_options:
    for feature in numerical_features:
        fig1 = plt.figure()
        sns.histplot(data=data, x=feature, color="black") 
        col1.pyplot(fig1)

# # if box plot selected on the multi select, creat box plot for each categorical features
if  'BOX-PLOT' in graph_options:
    for feature in categorical_features:
        fig2 = plt.figure()
        sns.boxplot( x=feature, y='PRICE',data=data, palette="Set1")
        # sns.catplot( data=data,x=feature, y='PRICE', kind="boxen")
        col2.pyplot(fig2)

# if count plot selected on the mulitselect, create count plot for each categorical features
if  'COUNT-PLOT' in graph_options:
    for feature in categorical_features:
        fig3 = plt.figure()
        sns.countplot(x=feature, data=data, palette="Set3")
        col3.pyplot(fig3)

