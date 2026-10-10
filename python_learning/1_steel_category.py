
import requests
from bs4 import BeautifulSoup
import joblib

response = requests.get("https://www.makeitfrom.com/material-group/Iron-Alloy")

bs = BeautifulSoup(response.text.encode(), features="html.parser")

steel_catos = bs.select("a[href*='Steel']")

joblib.dump(steel_catos, "steel_catos.jb")
