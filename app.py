from flask import Flask, render_template, request, flash, redirect, send_file, url_for, Request, send_from_directory
import os
from werkzeug.utils import secure_filename
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import pickle
import seaborn as sns
import subprocess
import ssri as ss
import eda_activity as eda_a
import pandas as pd

from pyfladesk import init_gui

UPLOAD_FOLDER = ''
ALLOWED_EXTENSIONS = {'csv', 'xls', 'xlsx'}


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 20 * 1000 * 1000
app.config['SECRET_KEY'] = '0a3eccba1fdadd919194172e'


def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route("/")
def hello():

    squery = ss.searchQuery("ssr")

    return render_template("index.html", squery=squery.head(10).to_html(classes="table table-striped"),
                           squery_tail=squery.tail(10).to_html(classes="table table-striped"), squery_shape=squery.shape)


@app.route("/", methods=["GET", "POST"])
def searchTarget():
    if request.method == "POST":
        squery = ss.searchQuery("ssr")
        if request.form.get('search') == "Search":
            nfrows = "First 10 rows"
            nlrows = "Last 10 rows"
            squery = ss.searchQuery("ssr")
            select_target = request.form.get('chem_id')
            starget = ss.selectTarget(select_target)
            pone = ss.phase_i()
            ptwo = eda_a.phaseii_Lipinski()
            return render_template("index.html", squery=squery.head(10).to_html(classes="table table-striped"), squery_tail=squery.tail(10).to_html(classes="table table-striped"), squery_shape=squery.shape,
                                   nfrows=nfrows, nlrows=nlrows, starget=starget.head(10).to_html(classes="table table-striped"), starget_tail=starget.tail(10).to_html(classes="table table-striped"), starget_shape=starget.shape,
                                   pone=pone.head(10).to_html(classes="table table-striped"),  pone_tail=pone.tail(10).to_html(classes="table table-striped"), pone_shape=pone.shape,
                                   ptwo=ptwo.head(10).to_html(classes="table table-striped"), ptwo_tail=ptwo.tail(10).to_html(classes="table table-striped"), ptwo_shape=ptwo.shape)

        elif request.form.get('upload_file') == "Upload":
            if 'file' not in request.files:
                flash('No file part')
                return redirect(request.url)
            file = request.files['file']
            # If the user does not select a file, the browser submits an
            # empty file without a filename.
            if file.filename == '':
                flash('No selected file')
                return redirect(request.url)
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                filename = "ssri_bioactivity.csv"
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                flash(filename + " has been uploaded successfully")
                return redirect('/downloadfile/' + filename)

            else:
                flash('Unsupported file format')
                return redirect(request.url)

        elif request.form.get('predict') == 'Predict':

            # Molecular descriptor calculator
            nfrows = "First 10 rows"
            nlrows = "Last 10 rows"
            df = pd.read_csv('ssri_bioactivity.csv')
            selection = ['canonical_smiles', 'mol_cid']
            df_selection = df[selection]
            df_selection.to_csv('df_selection.csv')
            df_selection.to_csv('molecule.smi', sep='\t',
                                index=False, header=False)

            def desc_calc():
                # Performs the descriptor calculation
                bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
                process = subprocess.Popen(
                    bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                os.remove('molecule.smi')

            # Molecular descriptor calculator
            # desc_calc()

            #   read descriptor out put produced by PADEL descriptor
            df_descriptor = pd.read_csv('descriptors_output.csv')

            # DETERMINE THE  X    V A R I A B L E from the descriptor output
            df_descriptor_x = df_descriptor.drop(columns=['Name'])

            # DETERMINE THE  Y    V A R I A B L E
            df_descriptor_y = df.iloc[:, -1]
            df_descriptor_y

            # combine X and Y variable
            dfCombined = pd.concat([df_descriptor_x, df_descriptor_y], axis=1)
            dfCombined.to_csv('descriptor_x_and_pIC50_y.csv', index=False)
            flash('descriptor_x_and_pIC50_y.csv saved to local storage')

            # remove low variance data to finally get the X independent variable
            def removeLowVariance(input_file, threshold=(.8 * (1 - .8))):
                selection = VarianceThreshold(threshold)
                selection.fit(input_file)
                return input_file[input_file.columns[selection.get_support(indices=True)]]

            Xvariance = removeLowVariance(
                df_descriptor_x, threshold=(.8 * (1 - .8)))
            Xaxis = pd.DataFrame(Xvariance)

            Xaxis.to_csv('descriptors_list.csv', index=False)
            flash('descriptors_list.csv saved to local storage')

            # finishedXaxis = Xaxis.fit_transform(df_descriptor)
            Xlist = list(pd.read_csv('descriptors_list.csv').columns)
            finishedXaxis = df_descriptor[Xlist]

            # Apply trained model to make prediction on query compounds

            def build_model(input_data):

                # Reads in saved regression model
                load_model = pickle.load(
                    open('ssri_bioactivity_model.pkl', 'rb'))
                # Apply model to make predictions
                prediction = load_model.predict(input_data)

                read_selection = pd.read_csv('df_selection.csv')
                prediction_output = pd.Series(prediction, name='pIC50')

                predictedVariables = pd.concat(
                    [read_selection, prediction_output], axis=1)
                # os.remove('predicted_outcome.csv')
                predictedVariables.to_csv("predicted_outcome.csv")
                return predictedVariables

            dof = build_model(finishedXaxis)
            return render_template("ssri.html", nfrows=nfrows, nlrows=nlrows, df=df.head(10).to_html(classes="table table-striped"), dfshape=df.shape, df_tail=df.tail(10).to_html(classes="table table-striped"),
                                   df_descriptor_x=df_descriptor_x.head(10).to_html(classes="table table-striped"), df_descriptor_x_tail=df_descriptor_x.tail(10).to_html(classes="table table-striped"), df_descriptor_xshape=df_descriptor_x.shape,
                                   df_selection=df_selection.head(10).to_html(classes="table table-striped"), df_selection_tail=df_selection.tail(10).to_html(classes="table table-striped"), df_selectionShape=df_selection.shape,
                                   df_descriptor=df_descriptor.head(10).to_html(classes="table table-striped"), df_descriptor_tail=df_descriptor.tail(10).to_html(classes="table table-striped"), df_des_shape=df_descriptor.shape,
                                   dfCombined=dfCombined.head(10).to_html(classes="table table-striped"), dfCombined_tail=dfCombined.tail(10).to_html(classes="table table-striped"), dfCombinedshape=dfCombined.shape,
                                   dof=dof.head(10).to_html(classes="table table-striped"), dof_tail=dof.tail(10).to_html(classes="table table-striped"), dofShape=dof.shape,
                                   X_variable=Xaxis.head(10).to_html(classes="table table-striped"),  X_variable_tail=Xaxis.tail(10).to_html(classes="table table-striped"), X_shape=Xaxis.shape)

    else:
        return render_template("index.html")


@app.route('/download-files/predicted_outcome.csv')
def return_files_tut():
    file_path = UPLOAD_FOLDER + 'predicted_outcome.csv'
    return send_file(file_path, as_attachment=True, attachment_filename='')


@app.route('/ssri')
def ssriHTML():
    return render_template('ssri.html')


if __name__ == '__main__':
    init_gui(app, window_title="SSRI MODEL")
    # app.run(port="200", debug=True)
