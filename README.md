# Procesamiento de secuencias y asignación taxonómica

Este repositorio contiene el pipeline utilizado para el procesamiento de secuencias y su posterior clasificación taxonómica.

## Preprocesamiento

- **Búsqueda y corte de primers**: realizado con [Cutadapt](https://github.com/marcelm/cutadapt).

- **Merge de secuencias**: realizado con [Flash2](https://github.com/dstreett/FLASH2).
  -- Min overlap 10
  -- Max overlap 100
    
- **Limpieza y clustering**:
  - Limpieza de secuencias con `UNOISE`.
  - Eliminación de *singletons* (parámetro --minsize 2 de UNOISE).
  - Agrupamiento de secuencias utilizando [SWARM](https://github.com/torognes/swarm).

## Clasificación taxonómica

Se consideraron varias estrategias de asignación:

- **Clasificación con Sintax**:
  - Usando la base de datos **Midori2** con un umbral de probabilidad `p > 0.7`.
  - Usando la base de datos **rCRUX**.

- **Clasificación con BLAST**:
  - Contra la base de datos **Midori2**.
  - Contra la base de datos **rCRUX**.

- **Consultas a BOLD**:
  - Para las secuencias que no obtienen asignación con Sintax o BLAST, se realizarán consultas a [BOLD Systems](https://boldsystems.org/) en grupos de 50 secuencias.

---

Este flujo de trabajo está orientado al análisis de datos de metabarcoding y permite comparar el rendimiento de diferentes enfoques de clasificación taxonómica.

