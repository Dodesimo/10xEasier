from flask import Flask, render_template, url_for, request
from tenxhelper import download

app = Flask(__name__)

@app.route("/")
def home():
    return render_template('index.html')

@app.route("/results",methods = ["POST","GET"])
def results(): 
    if request.method=="POST":
        download()
        return render_template('results.html')
