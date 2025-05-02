from django import forms


class genomeForm(forms.Form):
    genome_text = forms.CharField(widget=forms.Textarea(attrs={"rows":10,
                                                               "cols":60,
                                                               "placeholder": "MPXV genome(s) in .FASTA format"}
    ),
    required=False)

    uploaded_file = forms.FileField(required=False)