Welcome to the project to develop a computational system for immunological repertoires in case of celiac disease

In order to be able to download the project and make it work, please follow the following instructions (Further down we have added a youtube video to explain the steps clearly):

  1. Download it from github directly in zip format: 'ProjectLab-main.zip', for example in your Desktop.
  2. Extract all file from zip file in a new folder named 'ProjectLab-main
  3. Add the 'vdjbase_data' as well as the 'P1_CELIAC_METADATA.csv' file that you provided to us at the beginning of the final project in the folder where you extracted the codes from the zip file.
  4. Open the terminal in the folder where you extracted the project.
  5. Run the following instruction: 'pip3 install -r requirements.txt' in order to install all requires package for our project. Django is included inside this installation.
  6. Run instruction 'python3 manage.py makemigrations'.
  7. Run instruction 'python3 manage.py migrate'
  8. Run instruction 'python3 manage.py runserver'
  9. Open a new page on yout browser and enter the url: 'http://127.0.0.1:8000/data/home/'
  
 For now, you can see the home page of our interface, and then surf our interface. In our project book we explained how to use the interface and what each page consists of.
  
 In order to facilitate the installation of our project on your computer, here is a youtube video that we have prepared: https://www.youtube.com/watch?v=byq_7YUaB90
 
 Thank you for following the instructions. 
