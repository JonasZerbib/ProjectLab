from django.urls import path
from . import views

app_name = 'data'
urlpatterns = [
    path('home/', views.HomeView, name='home'),
    path('cdr3/', views.Cdr3View, name='cdr3'),
    path('cdr3results/', views.Cdr3Results, name='cdr3Results'),
    path('crossValidation/', views.CrossValidationView, name='crossValidation')
]
