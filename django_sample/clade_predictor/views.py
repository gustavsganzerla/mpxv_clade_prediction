from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.

def test_view(request):
    return HttpResponse('hi')

def home(request):
    return render(request, 'clade_predictor/home.html')