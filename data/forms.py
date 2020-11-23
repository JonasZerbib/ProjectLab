from django import forms
from django.forms import MultipleChoiceField

FEATURES_CHOICE = (('V_CALL', 'V_CALL'), ('D_CALL', 'D_CALL'), ('J_CALL', 'J_CALL'),
                   ('SEQUENCE_VDJ', 'SEQUENCE_VDJ'), ('SEQUENCE_IMGT', 'SEQUENCE_IMGT'),
                   ('V_SEQ_LENGTH', 'V_SEQ_LENGTH'),
                   ('JUNCTION_LENGTH', 'JUNCTION_LENGTH'), ('SEQUENCE_ID', 'SEQUENCE_ID'),
                   ('SEQUENCE_INPUT', 'SEQUENCE_INPUT'),
                   ('SEQUENCE_INPUT', 'SEQUENCE_INPUT'), ('FUNCTIONAL', 'FUNCTIONAL'), ('IN_FRAME', 'IN_FRAME'),
                   ('STOP', 'STOP'), ('MUTATED_INVARIANT', 'MUTATED_INVARIANT'), ('INDELS', 'INDELS'),
                   ('V_SEQ_START', 'V_SEQ_START'), ('V_GERM_START_VDJ', 'V_GERM_START_VDJ'),
                   ('V_GERM_LENGTH_VDJ', 'V_GERM_LENGTH_VDJ'),
                   ('V_GERM_START_IMGT', 'V_GERM_START_IMGT'), ('V_GERM_LENGTH_IMGT', 'V_GERM_LENGTH_IMGT'),
                   ('NP1_LENGTH', 'NP1_LENGTH'),
                   ('D_SEQ_START', 'D_SEQ_START'), ('D_SEQ_LENGTH', 'D_SEQ_LENGTH'), ('D_GERM_START', 'D_GERM_START'),
                   ('D_GERM_LENGTH', 'D_GERM_LENGTH'), ('NP2_LENGTH', 'NP2_LENGTH'), ('J_SEQ_START', 'J_SEQ_START'),
                   ('J_SEQ_LENGTH', 'J_SEQ_LENGTH'), ('J_GERM_START', 'J_GERM_START'),
                   ('J_GERM_LENGTH', 'J_GERM_LENGTH'),
                   ('JUNCTION', 'JUNCTION'), ('GERMLINE_IMGT', 'GERMLINE_IMGT'), ('FWR1_IMGT', 'FWR1_IMGT'),
                   ('FWR2_IMGT', 'FWR2_IMGT'), ('FWR3_IMGT', 'FWR3_IMGT'), ('FWR4_IMGT', 'FWR4_IMGT'),
                   ('CDR1_IMGT', 'CDR1_IMGT'), ('CDR2_IMGT', 'CDR2_IMGT'), ('CDR3_IMGT', 'CDR3_IMGT'),
                   ('CDR3_IGBLAST', 'CDR3_IGBLAST'), ('CDR3_IGBLAST_AA', 'CDR3_IGBLAST_AA'), ('PRCONS', 'PRCONS'),
                   ('SEQORIENT', 'SEQORIENT'), ('PRIMER', 'PRIMER'), ('ISOTYPE', 'ISOTYPE'),
                   ('CONSCOUNT', 'CONSCOUNT'), ('DUPCOUNT', 'DUPCOUNT'), ('DUPCOUNT_NEW', 'DUPCOUNT_NEW'),
                   ('CLONE', 'CLONE'), ('GERMLINE_V_CALL', 'GERMLINE_V_CALL'), ('GERMLINE_D_CALL', 'GERMLINE_D_CALL'),
                   ('GERMLINE_J_CALL', 'GERMLINE_J_CALL'), ('MUT', 'MUT'), ('CLONE_SIZE', 'CLONE_SIZE')
                   )

ALGO_CHOICE = (("Logistic regression", "Logistic regression"),
               ("Linear regression", "Linear regression"),
               ("SVM (Support Vector Machine)", "SVM (Support Vector Machine)"))

LABEL_CHOICE = (("age", "age"),
                ("helth statuts", "health statut"),
                ("sexe", "sex"))

DATASET = (("Celiac", "Celiac"),
           ("other", "other"))

CHECK_MODEL = (("Check model with One holdout cross validation", "Check model with One holdout cross validation"),
           ("Test your right model with machine learning algorithm", "Test your right model with machine learning "
                                                                     "algorithm"))

K_MERS_CHOICE = ((4, 4), (5, 5), (6, 6), (7, 7))

SPLIT_CHOICES = ((60, 60), (65, 65), (70, 70), (75, 75), (80, 80), (85, 85), (90, 90))


class CeliacForm(forms.Form):
    dataset = forms.ChoiceField(choices=DATASET, widget=forms.Select, label="")

class Cdr3Form(forms.Form):
    learning_rate = forms.FloatField(label="Enter your learning rate (preferably 0.01):")
    epochs = forms.IntegerField(label="Enter the number of epochs you want:")
    threshold = forms.FloatField(label="Enter the number of threshold you want:")
    seed = forms.IntegerField(label="Enter the number of seed you want:")
    kmers = forms.ChoiceField(choices=K_MERS_CHOICE, widget=forms.Select, label="Choose your kmers (lenght of snippet):")
    split = forms.ChoiceField(choices=SPLIT_CHOICES, label="How many percent of dataset you want for split")

    # check the value of learning rate
    def clean_learning_rate(self):
        learning_rate = self.cleaned_data['learning_rate']
        if learning_rate <= 0:
            raise forms.ValidationError("Please enter a positive number")
        return learning_rate

    # check the value of epochs
    def clean_epochs(self):
        epochs = self.cleaned_data['epochs']
        if epochs <= 0:
            raise forms.ValidationError("Please enter a positive number")
        return epochs

        # check the value threshold

    def clean_threshold(self):
        threshold = self.cleaned_data['threshold']
        if threshold <= 0 or threshold >= 1:
            raise forms.ValidationError("Please enter a number between 0 and 1")
        return threshold

class CrossValidationForm(forms.Form):
    learning_rate = forms.FloatField(label="Enter your learning rate (preferably 0.01):")
    epochs = forms.IntegerField(label="Enter the number of epochs you want:")
    threshold = forms.FloatField(label="Enter the number of threshold you want:")
    seed = forms.IntegerField(label="Enter the number of seed you want:")
    kmers = forms.ChoiceField(choices=K_MERS_CHOICE, widget=forms.Select, label="Choose your kmers (lenght of snippet):")

    # check_model = forms.ChoiceField(choices=CHECK_MODEL, widget=forms.Select, label="Test some model or test your right model if you already have the good one:")

    # check the value of learning rate
    def clean_learning_rate(self):
        learning_rate = self.cleaned_data['learning_rate']
        if learning_rate <= 0:
            raise forms.ValidationError("Please enter a positive number")
        return learning_rate

    # check the value of epochs
    def clean_epochs(self):
        epochs = self.cleaned_data['epochs']
        if epochs <= 0:
            raise forms.ValidationError("Please enter a positive number")
        return epochs

        # check the value threshold

    def clean_threshold(self):
        threshold = self.cleaned_data['threshold']
        if threshold <= 0 or threshold >= 1:
            raise forms.ValidationError("Please enter a number between 0 and 1")
        return threshold
