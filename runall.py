# Importações
from datetime import datetime
from time import time, sleep
from os import system
import re

# Parâmetros
n_runs = 3
command = './bios teste'
algorithms = [
   'StdPSO',
]
problems = [
   'TenBarTrussFrequencyFAST',
   'TenBarTrussFrequencyDIANA',
]

# Funções Utilitárias
def log(msg):
   r_file = open('results.txt', 'a')
   r_file.write(msg)
   r_file.close()

def replace(label, new_value):
   opt_file = open('teste.opt' , 'r')
   opt_data = opt_file.read()
   opt_file.close()
   opt_file = open('teste.opt' , 'w')
   opt_data = re.sub(f"%{label}\n\'.+\'", f"%{label}\n'{new_value}'", opt_data)
   opt_file.write(opt_data)
   opt_file.close()

def timed_log(msg):
   line = '-' * 50 + '\n'
   log(f'{line}{msg}\n{datetime.now()}\n{line}')

for a in algorithms:
   replace('OPTIMIZATION.ALGORITHM', a)
   for p in problems:
      replace('PROBLEM.TYPE', p)
      for i in range(n_runs):
         # Log do Teste
         timed_log(f'Algorithm: {a} | Problem: {p}\nRun: {i}')
      
         # Executando e Marcando Tempo
         t1 = time()
         system(command)
         t2 = time()
         log(f'Time Spent: {(t2 - t1):.2f} seconds\n')

         # Lendo Arquivo de Saída da Otimização
         out_file = open('teste.out', 'r')
         out_data = out_file.read()
         out_file.close()

         # Recuperando e Escrevendo Valores do Teste
         best_individuals = re.findall("%RESULT.BEST.INDIVIDUALS\n([^%]+)", out_data)[0]
         log(f'\nBest Individuals:\n{best_individuals.strip()}\n')
         f_obj = re.findall("%RESULT.OBJECTIVE.FUNCTION\n([^%]+)", out_data)[-1]
         log(f'\nObjective Function: {f_obj.strip()}\n')
         variables = re.findall("%DESIGN.VARIABLES\n([^%]+)", out_data)[-1]
         log(f'\nDesign Variables:\n{variables.strip()}\n')
         constraints = re.findall("%CONSTRAINT.VALUES\n([^%]+)", out_data)[-1]
         log(f'\nConstraints:\n{constraints.strip()}\n\n')

         # Esperando Próximo Minuto
         sleep(60)