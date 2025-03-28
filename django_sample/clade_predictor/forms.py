from django import forms


class genomeForm(forms.Form):
    genome_text = forms.CharField(widget=forms.Textarea(attrs={"rows":10,
                                                               "cols":60,
                                                               "placeholder": "MPXV genome(s) in .FASTA format"}
    ))

    uploaded_file = forms.FileField(
        label = 'Upload your genome(s) in .FASTA format',
        required=False
    )