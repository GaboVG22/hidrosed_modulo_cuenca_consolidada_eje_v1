# HidroSed · Cuenca Consolidada + Eje Activo v1

Módulo independiente para trabajar desde una **cuenca consolidada validada**.

## Qué hace

- Lee un KMZ/KML consolidado con cuenca de descarga, cuenca hidrológica, eje de cauce, puntos y curvas.
- Calcula **intercuenca = cuenca descarga - cuenca hidrológica**.
- Asocia el área incremental al tramo activo entre PC-HIDRO y PC-DESCARGA.
- Recorta curvas al corredor del eje.
- Genera curvas auxiliares interpoladas a lo largo del eje.
- Genera tabla de caudal incremental distribuido por km.
- Exporta KMZ, Excel, CSV y JSON.

## Ejecución

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Nota técnica

Las curvas auxiliares interpoladas son apoyo para preparación hidráulica. No son topografía levantada.
