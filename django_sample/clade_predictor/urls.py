from django.urls import path
from . import views

app_name = 'clade_predictor'

urlpatterns = [
    path("test_view/", views.test_view),
    path("home/", views.home, name = 'home')
    ]