from datetime import datetime
from siphon.simplewebservice.wyoming import WyomingUpperAir
from metpy.units import units, pandas_dataframe_to_unit_arrays
import matplotlib.pyplot as plt
import metpy.plots as plots

date = datetime(2023,9,10,0)
station = 'SLC'
df = WyomingUpperAir.request_data(date, station)
d = pandas_dataframe_to_unit_arrays(df) # turns to unidata arrays
print(d)
p = df['pressure'].values * units(df.units['pressure'])
T = df['temperature'].values * units(df.units['temperature'])
Td = df['dewpoint'].values * units(df.units['dewpoint'])
u= df['u_wind'].values * units(df.units['u_wind'])
v = df['v_wind'].values * units(df.units['v_wind'])

print(df.head())

fig = plt.figure(figsize=(10,10))
skew = plots.SkewT(fig)

skew.plot(p,T,'red')
skew.plot(p,Td,'green')
skew.plot_barbs(p,u,v)


fig.show()
skew.ax.set_ylim(150,20)
skew.ax.set_xlim(-80,0)

for p, t, h in zip(d['pressure'],d['temperature'],d['height']):
    if p <= 150 * units.hPa and p >= 20 * units.hPa:
        ft = h * units.feet * 3.2708
        skew.ax.text(1.06,p, round(ft.m,0), transform=skew.ax.get_yaxis_transform(which='tick2'))
skew.ax.text(1.06,157, '(feet)', transform=skew.ax.get_yaxis_transform(which='tick2'))
# set limits for pressure and temperature variables
fig.show()
xoffset = -.042
skew.ax.text(xoffset,150,'150',transform=skew.ax.get_yaxis_transform(which='tick2'))
skew.ax.text(xoffset,125,'125',transform=skew.ax.get_yaxis_transform(which='tick2'))
skew.ax.text(xoffset+0.012, 70,'70',transform=skew.ax.get_yaxis_transform(which='tick2'))
skew.ax.text(xoffset+0.012,50,'50',transform=skew.ax.get_yaxis_transform(which='tick2'))
skew.ax.text(xoffset+0.012,30,'30',transform=skew.ax.get_yaxis_transform(which='tick2'))
skew.ax.text(xoffset+0.012,20,'20',transform=skew.ax.get_yaxis_transform(which='tick2'))
skew.ax.axhline(y=125, xmin=0, xmax=1)
skew.ax.axhline(y=100, xmin=0, xmax=1)
skew.ax.axhline(y=70, xmin=0, xmax=1)
skew.ax.axhline(y=50, xmin=0, xmax=1)
skew.ax.axhline(y=30, xmin=0, xmax=1)
skew.ax.axhline(y=20, xmin=0, xmax=1)

fig.show()


print("finished base")

