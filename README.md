<a id="readme-top"></a>

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <h3 align="center">Geodesic on Sphere Visualization</h3>

  <p align="center">
    Geodesic on Sphere Visualization by using manim
    <br />
    <a href="https://github.com/hhaewon/sphere-geodesic/blob/main/Demo.mp4">View Demo</a>
    ·
    <a href="https://github.com/hhaewon/sphere-geodesic/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/hhaewon/sphere-geodesic/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
    </li>
    <li><a href="#usage">Usage</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

학교에서 평면 속에서 최단 거리는 직선이라고 배운다. 하지만 현실에서는 평면보다는 곡면이 더 많다. 그래서 곡면에서는 무엇이 최단거리인지 궁금했다. 가장 간단한 곡면 중 하나이며, 우리가 사는 지구의 모양이기도 한 구에서의 최단거리를 알아보고 싶었다. 

우리가 쓰는 세계지도는 메르카토르 도법에 의해 그려진다. 그러면 구의 두 지점 사이의 최단거리가 세계지도 상의 직선이라고 생각되기 쉽지만, 사실은 대원이 두 지점 사이의 최단거리이다. 이를 시각화 하기 위해, 세계지도 상의 직선을 구에 그리고 이를 대원과 비교했다. 또한 실제 길이도 나타내 정말 대원이 최단거리임을 숫자로 보여주었다.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Built With

* [manim][Manim-url]

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started
### Prerequisites

* [manim](https://docs.manim.community/en/stable/installation.html)

### Usage
change latitude and longitude in main.py
```python
### main.py
point1 = SphericalPoint(45, 90)
point2 = SphericalPoint(45, 0)
```
```sh
### shell
manim -pqh main.py SphereWithGeodesicScene 
```

### Execute in Colab
You can execute this in [this link](https://colab.research.google.com/drive/1bUQT0Bv52dVvW04c3DOCbvDKTM4Jjd3h?usp=sharing)  
Before you run all cells, you must upload [world_map.jpg](https://github.com/hhaewon/sphere-geodesic/blob/main/world_map.jpg) on colab

[Manim-url]: https://www.manim.community/
