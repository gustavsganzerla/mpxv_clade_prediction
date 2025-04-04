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
    model_clade = xgb.XGBClassifier()
    model_subclade = xgb.XGBClassifier()
    output = []
    n_seqs = 0

    model_clade.load_model('/var/www/django_app/django_sample/models/xgb_model.json')
    model_subclade.load_model('/var/www/django_app/django_sample/models/xgb_model_subclade.json')

    #model_clade.load_model('/Users/gustavosganzerla/mpxv_clade_prediction/django_sample/models/xgb_model.json')
    #model_subclade.load_model('/Users/gustavosganzerla/mpxv_clade_prediction/django_sample/models/xgb_model_subclade.json')


    
    if request.method == 'POST':
        form = genomeForm(request.POST)
        all_kmers = [''.join(kmer) for kmer in itertools.product('ATCG', repeat=3)]

        if form.is_valid():
            collected_data = form.cleaned_data
            genome_text = collected_data.get('genome_text')
            all_counts = []

            if genome_text:
               
                for record in SeqIO.parse(StringIO(genome_text), 'fasta'):
                   n_seqs +=1
                   kmer_counts = generate_kmers(record.seq, 3)
                   features = [kmer_counts.get(kmer, 0) for kmer in all_kmers]

                   X = np.array(features).reshape(1,-1)

                   
                   y_pred_proba = model_clade.predict_proba(X)
                   elements = y_pred_proba.tolist()
                   
                   for inner_list in elements:
                       clade = ''
                       if inner_list[0] > inner_list[1]:
                           clade = 'Clade 1'
                           y_pred_proba = model_subclade.predict_proba(X)
                           elements_subclade = y_pred_proba.tolist()

                           for inner_list_subclade in elements_subclade:
                               if inner_list_subclade[0] > inner_list_subclade[1]:
                                   clade += 'a'
                               else:
                                   clade += 'b'
                       else:
                           clade = 'Clade 2b'
                   output.append({
                       'id':record.id,
                       'prediction':clade,
                       'length':len(record.seq)
                   })
                    
                

                return render(request, 'clade_predictor/results.html', 
                              context={'output':output,
                                       'n_seqs':n_seqs})

    else:
        form = genomeForm()

        




    return render(request, 'clade_predictor/home.html',
                  context = {'form':form})