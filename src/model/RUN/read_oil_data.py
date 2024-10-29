import argparse
import pandas as pd

def read_oilbase(value,path,option):

   data = pd.read_csv(path,sep=';')

   if option=="NAME":
      try: #collecting line with correspondent oil name
         data.index = data.oil_name         
         oil = data.loc[value, :]

      except: #if there is no exact match
         # line below gets all names that start with the same letter as input
         oil_withletter = data[data['oil_name'].str.startswith(value[0])].oil_name.values.tolist()

         print('There is no oil with that name. Do you mean:')
         print('\n'.join(oil_withletter))
         
         raise ValueError('No oil with that name in database. Try suggestions mentioned above.')
   
   else: #uses closest API value to return an oil
      oil =  data.loc[(data.oil_api - float(value)).abs().idxmin()]

   #writing oil values to oil_file.txt that will be use in medslik-II simulation
   with open ('oil_file.txt','w') as f:
      f.write(f"{oil.oil_name}\n")
      f.write(f"{oil.oil_api:.02f}          API of Oil\n")
      f.write(f"{oil.oil_density:.03f}          Density of Oil\n")
      f.write(f"{oil.oil_res_dens:.03f}          Residual Density of Oil\n")
      f.write(f"{oil.oil_res_perc:.02f}          Residual Percent of Oil\n")
      f.write(f"{oil.oil_visc:.02f}          Viscosity of Oil\n")
      f.write(f"{oil.oil_temp:.02f}          Temperature at which Viscosity determined\n")
      f.write(f"{oil.oil_vap_press:.03f}          Vapour Pressure of Oil (bar)\n")
      
      f.close()

if __name__ == '__main__':

   parser = argparse.ArgumentParser(description='Collect oil information from the oilbase')
   parser.add_argument('option', type=str, help='Either "name" or "api"')
   parser.add_argument('value', help='Oil name or Oil API. Oil name must be exact value')
   

   args = parser.parse_args()

   path = 'oilbase.csv'

   read_oilbase(args.value,path,args.option)
