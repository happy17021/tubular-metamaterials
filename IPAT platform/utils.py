#-*- coding:utf-8 -*-
import glob
import re
import os.path
import sys

import cv2
import matplotlib.pyplot as plt
import numpy as np

def ShowImage(name, image):
    cv2.imshow(name, image)
    cv2.waitKey(0)  #
    cv2.destroyAllWindows()


def image2data(fname = './images/221.jpg',data_path = './data'):
    # fname = './images/221.jpg'
    image = cv2.imread(fname)
    grayIma = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)  #
    ret, image1 = cv2.threshold(grayIma, 127, 255, cv2.THRESH_BINARY)
    contours, hi = cv2.findContours(image1, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE) # opencv
    drawimage = image.copy()
    result1 = cv2.drawContours(drawimage, contours, -1, (0,255,0), 1)
    # ShowImage('result1', result1)
    #
    base = os.path.basename(fname)[:-4]
    c = 0
    for i,contour in enumerate(contours):
        if(cv2.contourArea(contour) <= 10):
            continue
        epsilon = 0.00000001 * cv2.arcLength(contour, True)
        approx = cv2.approxPolyDP(contour, epsilon, True)
        cv2.drawContours(image, [approx], 0, (255, 0, 0), 2)

        ap = [[i[0][0], i[0][1]] for i in approx]
        #
        scale = 12.8/128.
        ap = [(x * scale, y * scale) for (x, y) in ap]


        np.savetxt(os.path.join(data_path,f'{base}_{c+1}.txt'),ap)
        c += 1

    return os.path.basename(fname)[:-4],c


def update(fname = "Job-221.inp"):
    # fname = "Job-302.inp"
    with open(fname, 'r', encoding='utf-8', errors='ignore') as f:
        contents = f.read()
    new_content = ""
    f = True
    for i, line in enumerate(contents.split('\n')):
        if ("*Element, type=S3" in line):
            f = "NO"
            continue
        if ("*Element, type=S4R" in line):
            f = 'NO'
            continue
        if ("*Element, type=C3" in line):
            f = True
        if (f != "NO"):
            new_content += (line + "\n")
    with open(fname, 'w', encoding='utf-8', errors='ignore') as f:
        f.write(new_content)


def getAllImages(path = './images'):
    imgs = glob.glob(os.path.join(path,'*.jpg'))
    return imgs
