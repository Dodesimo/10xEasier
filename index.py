from flask import Flask, render_template, url_for, request
from tenxhelper import download, generateAnn, preprocess, qualityControl, generateGraphics, normalization, clustering

app = Flask(__name__)

@app.route("/")
def home():
    return render_template('index.html')

@app.route("/results",methods = ["POST","GET"])
def results(): 
    if request.method=="POST":
        download()
        data = generateAnn()
        data1 = preprocess(data=data)
        data2 = qualityControl(data=data1)
        generateGraphics(data=data2)
        data3 = normalization(data=data2)
        clustering(data=data3)

        return render_template('results.html')


