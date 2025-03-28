from django.urls import path
from . import views

urlpatterns = [
    path("test_view/", views.test_view)
    ]