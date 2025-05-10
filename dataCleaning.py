import pandas as pd
import os
print(os.getcwd())

dfBTC = pd.read_csv("data/BTCUSDT_1m.csv")
dfETH = pd.read_csv("data/ETHUSDT_1m.csv")
print("BTC Time Range:")
print("Start:", dfBTC["Open Time"].iloc[0])
print("End:  ", dfBTC["Open Time"].iloc[-1])
print(dfBTC.head())

print("\nETH Time Range:")
print("Start:", dfETH["Open Time"].iloc[0])
print("End:  ", dfETH["Open Time"].iloc[-1])
print(dfETH.head())
