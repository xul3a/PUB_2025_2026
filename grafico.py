import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "t10.tsv",
    sep="\t",
    comment="#",
    header=None,  # IMPORTANTE
    names=["t", "x", "y", "h", "K", "erro_rel"],
    index_col=False
)

#presas e predadores
plt.figure()
plt.plot(df["t"], df["x"], label="presas x(t)")
plt.plot(df["t"], df["y"], label="predadores y(t)")
plt.xlabel("tempo t")
plt.ylabel("população")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# plano de fase
plt.figure()
plt.plot(df["x"], df["y"])
plt.xlabel("presas x")
plt.ylabel("predadores y")
plt.title("Modelo Presa–Predador (Lotka–Volterra)")
plt.grid(True)
plt.tight_layout()
plt.show()

# erro relativo 
plt.figure()
plt.semilogy(df["t"], df["erro_rel"])
plt.xlabel("tempo t")
plt.ylabel("erro relativo de K")
plt.grid(True)
plt.tight_layout()
plt.show()

#passo adaptativo 
plt.figure()
plt.plot(df["t"], df["h"])
plt.xlabel("tempo t")
plt.ylabel("passo h")
plt.grid(True)
plt.tight_layout()
plt.show()
