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

from tensorflow.keras.models import load_model

# views to handle data
def generate_kmers(sequence, k):
    sequence = ''.join([base for base in sequence if base in 'ATCG'])

    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    return Counter(kmers)


def one_hot_encoding(seq):
    mapping = {'A':[1,0,0,0],
               'C':[0,1,0,0],
               'G':[0,0,1,0],
               'T':[0,0,0,1]}

    return np.array([mapping.get(base, [0,0,0,0]) for base in seq])

def split_genomes(sequence, window_size=1000, step_size=1000):
    windows = []
    for start in range(0, len(sequence) - window_size+1, step_size):
        window = sequence[start:start+window_size]
        windows.append(str(window))
    return windows


# Create your views here.
def home(request):
    model_clade = xgb.XGBClassifier()
    model_subclade = xgb.XGBClassifier()
    output = []
    n_seqs = 0

    model_clade.load_model('/var/www/django_app/django_sample/models/xgb_clade.json')
    model_subclade.load_model('/var/www/django_app/django_sample/models/xgb_model_subclade.json')

    model_partial_genomes_path = '/var/www/django_app/django_sample/models/v3-gpu_cnn_model.h5'
    #model_partial_genomes_path = '/Users/gustavosganzerla/mpxv_clade_prediction/django_sample/models/v3-gpu_cnn_model.h5'
    model_partial = load_model(model_partial_genomes_path)
    

    #model_clade.load_model('/Users/gustavosganzerla/mpxv_clade_prediction/django_sample/models/xgb_clade.json')
    #model_subclade.load_model('/Users/gustavosganzerla/mpxv_clade_prediction/django_sample/models/xgb_model_subclade.json')
    
    

    
    


    
    if request.method == 'POST':
        form = genomeForm(request.POST, request.FILES)
        all_kmers = [''.join(kmer) for kmer in itertools.product('ATCG', repeat=3)]

        if form.is_valid():
            collected_data = form.cleaned_data
            genome_text = collected_data.get('genome_text')
            uploaded_file = collected_data.get('uploaded_file')
            all_counts = []
            print(collected_data)


            ###here I am getting data from the text box in the website
            if genome_text:
                for record in SeqIO.parse(StringIO(genome_text), 'fasta'):
                    n_seqs +=1

                    ###classification with xgboost of complete genomes
                    if len(record.seq) >= 192000:
                        classification_type = 'Complete'
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
                    
                    elif len(record.seq) < 192000 and len(record.seq)>1000:
                        classification_type = 'Partial'
                        n_windows = 0
                        predictions_array = []
                        for window in split_genomes(str(record.seq)):
                            X_window = one_hot_encoding(window).flatten()
                            X_window = np.array(X_window)
                            X_window = X_window.reshape(-1, 1000, 4)
                            n_windows+=1
                            y_pred = model_partial.predict(X_window)
                            ###here, predictions_array will have a list of lists
                            predictions_array.append(np.argmax(y_pred, axis=1))
                            
                        ###here, predictions_array will have regular list
                        predictions_array = np.concatenate(predictions_array)
                        counts = np.bincount(predictions_array, minlength=3)
                        counts0, counts1, counts2 = counts[0], counts[1], counts[2]

                        if counts0 > counts1 and counts0 > counts2:
                            clade = 'Clade 1a'
                        elif counts1 > counts0 and counts1 > counts2:
                            clade = 'Clade 1b'
                        elif counts2 > counts0 and counts2 > counts1:
                            clade = 'Clade 2b'
                        
                    elif len(record.seq) < 1000:
                        classification_type = 'Input too short'
                        clade = 'Underdetermined'
                    
                    output.append({
                            'id':record.id,
                            'prediction':clade,
                            'length':len(record.seq),
                            'classification_type':classification_type
                        })
                    

                    

                return render(request, 'clade_predictor/results.html', 
                              context={'output':output,
                                       'n_seqs':n_seqs})

            
            if uploaded_file:
                uploaded_file_data = uploaded_file.read().decode('utf-8')
                uploaded_file_io = StringIO(uploaded_file_data)

                for record in SeqIO.parse(uploaded_file_io, 'fasta'):
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