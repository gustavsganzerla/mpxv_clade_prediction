from django.urls import path
from . import views

app_name = 'clade_predictor'

urlpatterns = [
    path("home/", views.home, name = 'home')
    ]