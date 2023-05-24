from django.shortcuts import render
from django.shortcuts import redirect
from django.http import HttpResponse
import os
import glob
import datetime
# Create your views here.
def redirect_view(request):
    response = redirect('/home')
    return response
def data(request):
    data_dir = '/home/evanz/data/sif_test/'
    files = glob.glob('%s/chunk*/' % data_dir)
    dates = [datetime.datetime.fromtimestamp(os.path.getmtime(file)) for file in files]
    num = list(range(1, len(files)+1))
    merge = tuple(zip(num, dates))
    context = {'files': merge}   
    return render(request, 'data.html', context=context)

def scan(request, scan_num):
    scan_dir = '/home/evanz/data/sif_test/chunk%d' % scan_num
    files = glob.glob('%s/*png' % scan_dir)
    files.sort(key=os.path.getmtime)  
    dates = [datetime.datetime.fromtimestamp(os.path.getmtime(file)) for file in files]
    files = ['chunk' + str(scan_num) + '/' + os.path.basename(x) for x in files]
    num = list(range(1, len(files)+1))
    merge = tuple(zip(num, files, dates))
    context = {'files': merge, 'scan_num': scan_num}
    return render(request, 'scan.html', context=context)

def image(request, scan_num, image_num):
    img_dir = '/home/evanz/data/sif_test/chunk%d' % scan_num
    imgs = glob.glob('%s/*png' % img_dir)
    imgs.sort(key=os.path.getmtime)
    imgs = ['chunk' + str(scan_num) + '/' + os.path.basename(x) for x in imgs]
    img_src = imgs[image_num-1]
    context = {'img_src': img_src}
    return render(request, 'image.html', context=context)
