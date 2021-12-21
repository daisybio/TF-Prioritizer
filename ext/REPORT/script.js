function accordion() {
    let acc = document.getElementsByClassName("accordion");
    let i;

    for (i = 0; i < acc.length; i++) {
        acc[i].addEventListener("click", function () {
            this.classList.toggle("active");

            let panel = this.nextElementSibling;
            if (panel.style.display === "grid") {
                panel.style.display = "none";
            } else {
                panel.style.display = "grid";
            }
        });
    }
}

function selector() {
    let selectors = document.getElementsByClassName("selector");
    let i;

    for (i = 0; i < selectors.length; i++) {
        selectors[i].addEventListener("click", function () {
            let setClasses = !this.classList.contains("active");
            setClass(selectors, "active", "remove");

            let image = document.getElementById(this.name)

            image.src = this.value + ".png";

            if (setClasses) {
                this.classList.toggle("active");
            }
        });
    }
}

function setClass(els, className, fnName) {
    for (let i = 0; i < els.length; i++) {
        els[i].classList[fnName](className);
    }
}