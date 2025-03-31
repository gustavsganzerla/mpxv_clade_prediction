from django.shortcuts import render
from django.http import HttpResponse
from . forms import genomeForm
from io import StringIO
from Bio import SeqIO
from collections import Counter
import itertools
import pandas as pd
import numpy as np
import xgboost as xgb

# views to handle data
def generate_kmers(sequence, k):
    sequence = ''.join([base for base in sequence if base in 'ATCG'])

    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    return Counter(kmers)



# Create your views here.
def home(request):
    loaded_model = xgb.XGBClassifier()
    output = []

    loaded_model.load_model('/Users/gustavosganzerla/mpxv_clade_prediction/django_sample/models/xgb_model.json')
    if request.method == 'POST':
        form = genomeForm(request.POST)
        all_kmers = [''.join(kmer) for kmer in itertools.product('ATCG', repeat=3)]

        if form.is_valid():
            collected_data = form.cleaned_data
            genome_text = collected_data.get('genome_text')
            all_counts = []

            if genome_text:
               
                for record in SeqIO.parse(StringIO(genome_text), 'fasta'):
                   kmer_counts = generate_kmers(record.seq, 3)
                   features = [kmer_counts.get(kmer, 0) for kmer in all_kmers]

                   X = np.array(features).reshape(1,-1)
                   y_pred_proba = loaded_model.predict_proba(X)
                   elements = y_pred_proba.tolist()
                   
                   for inner_list in elements:
                       clade = ''
                       if inner_list[0] > inner_list[1]:
                           clade = '1'
                       else:
                           clade = '2'
                   output.append({
                       'ID':record.id,
                       'Prediction':clade
                   })
                    
                print(output)

                return render(request, 'clade_predictor/results.html', 
                              context={'output':output})

    else:
        form = genomeForm()

        




    return render(request, 'clade_predictor/home.html',
                  context = {'form':form})