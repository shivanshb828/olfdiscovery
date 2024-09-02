import csv
import time
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Initialize the WebDriver using webdriver_manager
driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()))

# Open the MolPort website
driver.get("https://www.molport.com/shop/index")

def search_zinc_id(zinc_id):
    try:
        # Find the search box
        search_box = driver.find_element(By.ID, 'search-form__top-header_value')  # Adjust the search box ID if needed

        # Clear any existing text in the search box
        search_box.clear()

        # Enter the ZINC ID and submit the form
        search_box.send_keys(zinc_id)
        search_box.send_keys(Keys.RETURN)

        # Wait for the results page to load
        WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.CLASS_NAME, 'container'))  # Adjust as needed
        )

        # Check the current URL for the "no-matching-structures-were-found" string
        current_url = driver.current_url
        if "no-matching-structures-were-found" not in current_url:
            return True
        else:
            return False

    except Exception as e:
        print(f"Error searching for {zinc_id}: {e}")
        return False

# Open the CSV file and iterate through the ZINC IDs
with open('maestroResults.csv', mode='r') as file:
    csv_reader = csv.DictReader(file)
    with open('results.csv', mode='w', newline='') as results_file:
        fieldnames = ['ZINC ID']
        writer = csv.DictWriter(results_file, fieldnames=fieldnames)
        writer.writeheader()

        for row in csv_reader:
            zinc_id = row['title']
            print(f"Searching for {zinc_id}...")
            if search_zinc_id(zinc_id):
                writer.writerow({'ZINC ID': zinc_id})
                print(f"Added {zinc_id} to results.csv")


# Close the WebDriver
driver.quit()