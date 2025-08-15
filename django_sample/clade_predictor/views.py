from django.shortcuts import render
from django.http import HttpResponse
from . forms import genomeForm, contactForm
from io import StringIO
from Bio import SeqIO
from collections import Counter
import itertools
import pandas as pd
import numpy as np
import xgboost as xgb
from itertools import product

from tensorflow.keras.models import load_model


from django.core.mail import EmailMessage, get_connection
from django.conf import settings

# views to handle data
bases = ['A', 'C', 'G', 'T']
trimer_list = [''.join(p) for p in product(bases, repeat=3)]
WINDOW_SIZE = 1000
STEP_SIZE = 500
NUC_DICT = {'A': [1,0,0,0],
            'T': [0,1,0,0],
            'C': [0,0,1,0],
            'G': [0,0,0,1]}
###functions to process input data

def trinuc_freq(seq):
    seq = seq.upper()
    counts = {tri: 0 for tri in trimer_list}
    total = len(seq) - 2
    for i in range(total):
        tri = seq[i:i+3]
        if all(b in bases for b in tri):
            counts[tri] += 1
    freqs = np.array([counts[tri]/total if total > 0 else 0 for tri in trimer_list])
    return freqs


def one_hot_encode(seq):
    return np.array([NUC_DICT.get(base.upper(), [0,0,0,0]) for base in seq])



# Create your views here.
def home(request):
    model_clade = xgb.XGBClassifier()
    model_subclade = xgb.XGBClassifier()
    output = []
    n_seqs = 0

    ###MODELS
    ###COMPLETE GENOME
    MODEL_COMPLETE_CLADE_PATH = '/var/www/django_app/django_sample/models/complete_clade.json'
    MODEL_COMPLETE_SUBCLADE_PATH = '/var/www/django_app/django_sample/models/complete_subclade.json'

    model_complete_clade = xgb.XGBClassifier()
    model_complete_subclade = xgb.XGBClassifier()

    model_complete_clade.load_model(MODEL_COMPLETE_CLADE_PATH)
    model_complete_subclade.load_model(MODEL_COMPLETE_SUBCLADE_PATH)

    ###PARTIAL GENOME
    MODEL_PARTIAL_CLADE_PATH = '/var/www/django_app/django_sample/models/partial_clade.h5'
    MODEL_PARTIAL_SUBCLADE_PATH = '/var/www/django_app/django_sample/models/partial_subclade.h5'

    model_partial_clade = load_model(MODEL_PARTIAL_CLADE_PATH)
    model_partial_subclade = load_model(MODEL_PARTIAL_SUBCLADE_PATH)
    


    
    
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
                    seq = record.seq
                    n_seqs +=1
                    clade = ''

                    ###classification with xgboost of complete genomes
                    if len(seq) >= 188000:
                        classification_type = 'Complete (XGBoost)'
                        
                        freq_vector = trinuc_freq(seq)
                        df_features = pd.DataFrame([freq_vector], columns=trimer_list)

                        pred_clade = model_complete_clade.predict(df_features)[0]

                        if pred_clade == 1:
                            clade = 'Subclade IIb'
                        else:
                            pred_subclade = model_complete_subclade.predict(df_features)[0]

                            if pred_subclade == 0:
                                clade = 'Subclade Ia'
                            else:
                                clade = 'Subclade Ib'
                    
                    elif len(record.seq) < 192000 and len(record.seq)>1000:
                        classification_type = 'Partial (CNN)'
                        X = []

                        for start in range(0, len(seq) - WINDOW_SIZE + 1, STEP_SIZE):
                            window_seq = seq[start:start+WINDOW_SIZE]

                            if len(seq) == WINDOW_SIZE:
                                X.append(one_hot_encode(window_seq))
                        if X:
                            X = np.array(X)

                            pred_clade_array = model_partial_clade.predict(X)
                            pred_clade_array = np.argmax(pred_clade_array, axis=1)

                            clade1_votes = np.sum(pred_clade_array==0)
                            clade2_votes = np.sum(pred_clade_array==1)

                            if clade2_votes == clade1_votes:
                                clade = 'Undetermined'
                            if  clade2_votes > clade1_votes:
                                clade = 'Subclade IIb'

                            else:
                                pred_subclade_array = model_partial_subclade.predict(X)
                                pred_subclade_array = np.argmax(pred_subclade_array, axis=1)

                                clade1a_votes = np.sum(pred_subclade_array==0)
                                clade1b_votes = np.sum(pred_subclade_array==1)

                                if clade1a_votes == clade1b_votes:
                                    clade = 'Undetermined'
                                if clade1a_votes > clade1b_votes:
                                    clade = 'Subclade Ia'
                                else:
                                    clade = 'Subclade Ib'
                                

                        
                        
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

            ###here I am getting data from uploaded file in the website
            if uploaded_file:
                uploaded_file_data = uploaded_file.read().decode('utf-8')
                uploaded_file_io = StringIO(uploaded_file_data)

                for record in SeqIO.parse(StringIO(genome_text), 'fasta'):
                    seq = record.seq
                    n_seqs +=1
                    clade = ''

                    ###classification with xgboost of complete genomes
                    if len(seq) >= 188000:
                        classification_type = 'Complete (XGBoost)'
                        
                        freq_vector = trinuc_freq(seq)
                        df_features = pd.DataFrame([freq_vector], columns=trimer_list)

                        pred_clade = model_complete_clade.predict(df_features)[0]

                        if pred_clade == 1:
                            clade = 'Subclade IIb'
                        else:
                            pred_subclade = model_complete_subclade.predict(df_features)[0]

                            if pred_subclade == 0:
                                clade = 'Subclade Ia'
                            else:
                                clade = 'Subclade Ib'
                    
                    elif len(record.seq) < 192000 and len(record.seq)>1000:
                        classification_type = 'Partial (CNN)'
                        X = []

                        for start in range(0, len(seq) - WINDOW_SIZE + 1, STEP_SIZE):
                            window_seq = seq[start:start+WINDOW_SIZE]

                            if len(seq) == WINDOW_SIZE:
                                X.append(one_hot_encode(window_seq))
                        if X:
                            X = np.array(X)

                            pred_clade_array = model_partial_clade.predict(X)
                            pred_clade_array = np.argmax(pred_clade_array, axis=1)

                            clade1_votes = np.sum(pred_clade_array==0)
                            clade2_votes = np.sum(pred_clade_array==1)

                            if clade2_votes == clade1_votes:
                                clade = 'Undetermined'
                            if  clade2_votes > clade1_votes:
                                clade = 'Subclade IIb'

                            else:
                                pred_subclade_array = model_partial_subclade.predict(X)
                                pred_subclade_array = np.argmax(pred_subclade_array, axis=1)

                                clade1a_votes = np.sum(pred_subclade_array==0)
                                clade1b_votes = np.sum(pred_subclade_array==1)

                                if clade1a_votes == clade1b_votes:
                                    clade = 'Undetermined'
                                if clade1a_votes > clade1b_votes:
                                    clade = 'Subclade Ia'
                                else:
                                    clade = 'Subclade Ib'
                                

                        
                        
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
    else:
        form = genomeForm()



    return render(request, 'clade_predictor/home.html',
                  context = {'form':form})

def contact(request):
    if request.method=='POST':
        form = contactForm(request.POST)
        if form.is_valid():
            collected_data = form.cleaned_data
            subject = collected_data.get('subject')
            email = collected_data.get('email')
            message = collected_data.get('message')

            with get_connection(
                host = settings.EMAIL_HOST,
                port = settings.EMAIL_PORT,
                username = settings.EMAIL_HOST_USER,
                password = settings.EMAIL_HOST_PASSWORD,
                use_ssl = settings.EMAIL_USE_SSL
            ) as connection:
                subject = "MPP_"+subject
                email_from = settings.EMAIL_HOST_USER
                recipient_list = ['sganzerlagustavo@gmail.com']
                message = f"{message}\n{email}"

                email = EmailMessage(
                    subject,
                    message,
                    email_from,
                    recipient_list
                )
                email.send()
                return render(request, 'clade_predictor/contact_success.html')

    else:
        form = contactForm()

    
    return render(request, 'clade_predictor/contact.html',
                  {'form':form})
    
def how_to_use(request):
    return render(request, 'clade_predictor/how_to_use.html')