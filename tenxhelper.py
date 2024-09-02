from flask import request
import os
#def createAnnData(filename):
    

def download():
    downloads = request.files.getlist("uploadbutton")
    
    for download in downloads:
        download.save(os.path.join("/Users/atharvainamdar/Desktop/Appseq/App/downloads", download.filename))