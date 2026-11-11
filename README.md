# Reconstrucción Arqueológica mediante Spline Cúbico Paramétrico Natural

![HTML5](https://img.shields.io/badge/HTML5-E34F26?style=for-the-badge&logo=html5&logoColor=white)
![JavaScript](https://img.shields.io/badge/JavaScript-F7DF1E?style=for-the-badge&logo=javascript&logoColor=black)
![CSS3](https://img.shields.io/badge/CSS3-1572B6?style=for-the-badge&logo=css3&logoColor=white)
![Licencia](https://img.shields.io/badge/Licencia-MIT-green.style=for-the-badge)

Aplicación web interactiva desarrollada para la asignación de **Métodos Numéricos**, diseñada para la interpolación y reconstrucción visual del perfil de piezas arqueológicas sobre fotografías digitales. Utiliza **Splines Cúbicos Paramétricos Naturales** parametrizados por longitud de cuerda y resueltos mediante el **Algoritmo de Thomas**.

---

## Descripción del Proyecto

La herramienta permite cargar la fotografía de una pieza o fragmento arqueológico, colocar puntos de control sobre su contorno y calcular una curva suavizada e interpolada de alta precisión. Todo el cálculo matemático se realiza directamente en el navegador de manera autónoma y sin dependencias externas.

---

## Características Principales

* **Lienzo Interactivo (Canvas 2D):**
  * Carga de imágenes locales mediante botón o arrastrar y soltar (*drag & drop*).
  * Control interactivo de puntos: agregar, mover, eliminar uno a uno o limpiar lienzo.
  * Ajuste fino de coordenadas con teclado y zoom/desplazamiento en tiempo real.
* **Modos de Visualización:**
  * Renderizado de líneas auxiliares, coloreado por segmentos e inspección de coordenadas.
  * Rellenado del perfil reconstruido.
  * Comparador interactivo "antes / después" mediante un divisor deslizable.
* **Animación de Reconstrucción:**
  * Control de reproducción, pausa, ajuste de velocidad ($0.25\times$ a $3.00\times$) y barra de reproducción por tiempo.
* **Panel de Verificación Numérica:**
  * Inspección en tiempo real de los coeficientes polinómicos ($a, b, c, d$) para $X(t)$ e $Y(t)$ en cada intervalo.
  * Resumen estadístico: puntos de control, número de intervalos y longitud de cuerda total.
* **Exportación de Resultados:**
  * Exportación del perfil renderizado en formato **PNG**.
  * Exportación de los datos del modelo y puntos en formato **JSON**.
* **Funcionamiento Offline:** Desarrollado 100% en *Vanilla JavaScript*, HTML5 y CSS3. Funciona sin necesidad de servidor ni conexión a Internet.

---

## Fundamentos Matemáticos

Para evitar problemas en contornos cerrados o de múltiples valores donde una función simple $y = f(x)$ falla, la curva se parametriza respecto a la **longitud de cuerda acumulada** $t$.

### 1. Parametrización por Longitud de Cuerda
Dados $n+1$ puntos de control $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$:

$$t_0 = 0$$

$$t_i = t_{i-1} + \sqrt{(x_i - x_{i-1})^2 + (y_i - y_{i-1})^2}, \quad i = 1, \dots, n$$

$$h_i = t_{i+1} - t_i$$

---

### 2. Ecuación del Spline Cúbico por Tramo
Para cada intervalo $t \in [t_i, t_{i+1}]$ (para $i = 0, \dots, n-1$), las coordenadas $X(t)$ e $Y(t)$ se definen como:

$$X_i(t) = a_{x,i} + b_{x,i}(t - t_i) + c_{x,i}(t - t_i)^2 + d_{x,i}(t - t_i)^3$$

$$Y_i(t) = a_{y,i} + b_{y,i}(t - t_i) + c_{y,i}(t - t_i)^2 + d_{y,i}(t - t_i)^3$$

donde $a_{x,i} = x_i$ y $a_{y,i} = y_i$.

---

### 3. Condiciones de Continuidad y Frontera Natural

El algoritmo garantiza suave continuidad hasta la segunda derivada:

* **Continuidad $C^0$:** $X_i(t_{i+1}) = X_{i+1}(t_{i+1})$
* **Continuidad $C^1$:** $b_i + 2c_i h_i + 3d_i h_i^2 = b_{i+1}$
* **Continuidad $C^2$:** $2c_i + 6d_i h_i = 2c_{i+1}$
* **Condición de Frontera Natural:** $X''(t_0) = X''(t_n) = 0 \implies c_0 = c_n = 0$

El sistema de ecuaciones tridiagonal resultante para $c_i$:

$$h_{i-1}c_{i-1} + 2(h_{i-1} + h_i)c_i + h_i c_{i+1} = 3 \left( \frac{a_{i+1} - a_i}{h_i} - \frac{a_i - a_{i-1}}{h_{i-1}} \right)$$

se resuelve de forma eficiente mediante el **Algoritmo de Thomas** de manera independiente para la componente $X(t)$ y la componente $Y(t)$.

Posteriormente, los coeficientes restantes se calculan mediante:

$$b_i = \frac{a_{i+1} - a_i}{h_i} - \frac{h_i(2c_i + c_{i+1})}{3}$$

$$d_i = \frac{c_{i+1} - c_i}{3h_i}$$

---

## Estructura del Repositorio

```text
.
├── index.html       # Estructura principal de la aplicación e interfaz de usuario
├── estilos.css      # Hojas de estilo para la maquetación y controles
└── matematica.js    # Lógica de interacción, resolución tridiagonal y renderizado 2D
```

---

## Uso Rápido

1. Clona el repositorio o descarga los archivos:
   ```bash
   git clone https://github.com/tu-usuario/tu-repositorio.git
   ```
2. Abre el archivo `index.html` en cualquier navegador web moderno (Chrome, Firefox, Edge, Safari).
3. Haz clic en **"Cargar fotografía"** para subir la imagen de una pieza arqueológica.
4. Haz clic a lo largo del contorno de la pieza para definir los puntos de control (se recomiendan entre 8 y 15 puntos).
5. Presiona **"CALCULAR SPLINE"** para generar la reconstrucción y examinar las métricas numéricas.

---

## Controles e Interacción

| Acción | Control |
| :--- | :--- |
| **Agregar punto** | Clic izquierdo sobre la imagen |
| **Mover punto** | Arrastrar punto existente |
| **Eliminar punto** | Clic derecho sobre el punto |
| **Ajuste fino** | Teclas de dirección (Flechas). Mantener `Shift` para avanzar 5 px |
| **Zoom** | Rueda del ratón |
| **Desplazar lienzo** | `Shift` + arrastrar, clic con botón central o mantener tecla `Espacio` |

---

## Licencia

Este proyecto está disponible bajo la licencia **MIT**. Consulta el archivo `LICENSE` para más detalles.
