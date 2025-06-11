from django.urls import path
from . import views

app_name = 'clade_predictor'

urlpatterns = [
    path("home/", views.home, name = 'home'),
    path("contact/", views.contact, name = 'contact'),
    path("how_to_use/", views.how_to_use, name = 'how_to_use')
    ]