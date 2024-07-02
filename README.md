# bionetquery
BioNetQuery: Cell Type-Specific Biomarker and Pathway Analysis

## Commands
New virtual env:
```sh
virtualenv -p python3 .venv
```

Activate virtual env:
```sh
source .venv/bin/activatie
```

Install packages: 
```sh
pip install -r requirements.txt
```

See all installed packages:
```sh
pip freeze
```

Write out all installed packages to requirements.txt
```sh
pip freeze > requirements.txt
```

Run 'function one'
```sh
python src/python_functions/function_one.py podocyte
```