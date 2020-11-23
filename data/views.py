from random import *

from django.forms import forms
from django.shortcuts import render, redirect
from data import forms
from datasetCeliac.loadData import get_dataset_celiac
from datasetCeliac.modelLogisticRegression import train, test, initialization_of_parameters, split_data, crossVal
from datasetCeliac.presentation_of_results import presentation_of_results, cross_val_results
from datasetCeliac.presentationOfDataset import display_metaData
import numpy as np


# Create your views here.

def HomeView(request):
    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = forms.CeliacForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            dataset = form.cleaned_data["dataset"]
            request.session['dataset'] = dataset  # set 'dataset' in the session
            return redirect('data:cdr3')
        return redirect('home')
    else:  # if this is a GET request or the first time we access this view
        context = {'form': forms.CeliacForm(), "resumeOfDatasets": display_metaData()}
        return render(request, 'data/home.html', context)


def Cdr3View(request):
    if request.method == 'POST':  # if this is a POST request we need to process the form data
        # create a form instance and populate it with data from the request:
        my_form = forms.Cdr3Form(request.POST)
        # check whether it's valid:
        if my_form.is_valid():
            # process the data in form.cleaned_data as required
            kmers = int(my_form.cleaned_data['kmers'])  # recover kmers from user
            threshold = my_form.cleaned_data['threshold']  # recover threshold from user
            learning_rate = my_form.cleaned_data['learning_rate']  # recover learning rate from user
            epochs = my_form.cleaned_data['epochs']  # recover epochs from user
            split = my_form.cleaned_data['split']  # recover split from user
            seed = my_form.cleaned_data['seed']  # recover seed from user
            x_dataset, y_dataset = get_dataset_celiac(int(kmers))  # load dataset from datasetCeliac.loadData file
            train_x, train_y, test_x, test_y = split_data(int(split), x_dataset, y_dataset)  # split data
            parameters = initialization_of_parameters(int(kmers), int(seed))  # initialize parameters w, b0 and bfreq
            train_results = train(parameters, epochs, learning_rate, threshold, train_x, train_y)  # train the model
            test_results = test(test_x, test_y, threshold, train_results)  # test the model
            img1, img2, img3, img4, img6 = presentation_of_results(
                test_results)  # recover plots of results in order to display them in the next webpage
            return render(request, 'data/cdr3Results.html',
                          {'data1': img1, 'data2': img2, 'data3': img3, 'data4': img4, 'data5': img6})

    else:  # if a GET (or any other method) we'll create a blank form
        dataset = request.session['dataset']  # get 'dataset' from the session
        context = {'form': forms.Cdr3Form(), "dataset": dataset}
        return render(request, 'data/cdr3.html', context)


def CrossValidationView(request):
    # if this is a POST request we need to process the form data
    if request.method == 'POST':  # if this is a POST request we need to process the form data
        # create a form instance and populate it with data from the request:
        my_form = forms.CrossValidationForm(request.POST)
        # check whether it's valid:
        if my_form.is_valid():
            # process the data in form.cleaned_data as required
            kmers = int(my_form.cleaned_data['kmers'])  # recover kmers from user
            threshold = my_form.cleaned_data['threshold']  # recover threshold from user
            learning_rate = my_form.cleaned_data['learning_rate']  # recover learning rate from user
            epochs = my_form.cleaned_data['epochs']  # recover epochs from user
            seed = my_form.cleaned_data['seed']
            x_dataset, y_dataset = get_dataset_celiac(int(kmers))  # load dataset from datasetCeliac.loadData file
            parameters = initialization_of_parameters(int(kmers), int(seed))  # initialize parameters w, b0 and bfreq
            prediction, last_acc = crossVal(parameters, epochs, learning_rate, threshold, x_dataset,
                                            y_dataset)  # train the model
            img1 = cross_val_results(last_acc, prediction)  # recover plot of results in order to display them in the
            # next webpage
            return render(request, 'data/crossValidationResults.html', {'data1': img1})
    # if a GET (or any other method) we'll create a blank form
    else:
        dataset = request.session['dataset']  # get 'student_id' from the session
        context = {'form': forms.CrossValidationForm(), "dataset": dataset}
        return render(request, 'data/crossValidation.html', context)


def Cdr3Results(request):
    context = {'form': forms.Cdr3Form()}
    return render(request, 'data/cdr3.html', context)
