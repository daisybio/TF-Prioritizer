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

function openImageInTab(id) {
    let image = document.getElementById(id);
    window.open(image.src);
}

function openModal(modal_id) {
    let modal = document.getElementById(modal_id);
    modal.style.display = "block";
}

function closeModal(modal_id) {
    let modal = document.getElementById(modal_id);
    modal.style.display = "none";
}

function setClass(els, className, fnName) {
    for (let i = 0; i < els.length; i++) {
        els[i].classList[fnName](className);
    }
}

function tableMouseOver(gene, col, row) {
    colHead = document.getElementById(gene + "-col-" + col);
    rowHead = document.getElementById(gene + "-row-" + row);

    colHead.style.backgroundColor = "rgba(var(--black), 0.1)";
    rowHead.style.backgroundColor = "rgba(var(--black), 0.1)";
}

function tableMouseOut(gene, col, row) {
    colHead = document.getElementById(gene + "-col-" + col);
    rowHead = document.getElementById(gene + "-row-" + row);

    colHead.style.backgroundColor = "initial";
    rowHead.style.backgroundColor = "initial";
}

function init() {
    let hm_selectors = document.getElementsByClassName("hm-selector");
    hm_selectors[0].click();

    select_first_possible_group();

    update_image();
}

function select_first_possible_group() {
    let group_selectors = document.getElementsByClassName("group-selector");
    let i;

    for (i = 0; i < group_selectors.length; i++) {
        if (!group_selectors[i].disabled) {
            select_group(group_selectors[i]);
            break;
        }
    }
}

function hm_selector() {
    let hm_selectors = document.getElementsByClassName("hm-selector");
    let i;

    for (i = 0; i < hm_selectors.length; i++) {
        hm_selectors[i].addEventListener("click", function () {
            if (this.classList.contains("active")) return;

            let setClasses = !this.classList.contains("active");
            setClass(hm_selectors, "active", "remove");

            if (setClasses) {
                this.classList.toggle("active");
            }

            let possibleGroups = Object.keys(combinations[this.value]);

            let group_selectors = document.getElementsByClassName("group-selector");
            let j;

            for (j = 0; j < group_selectors.length; j++) {
                group_selectors[j].disabled = !possibleGroups.includes(group_selectors[j].value);
            }

            let activeGroup = get_active_element("group-selector");

            if (!(typeof activeGroup === 'undefined')) {
                if (!possibleGroups.includes(activeGroup.value)) {
                    select_first_possible_group();
                } else {
                    select_group(activeGroup);
                }
            }
        });
    }
}

function select_group(element) {
    let group_selectors = document.getElementsByClassName("group-selector");

    if (!element.classList.contains("active")) {
        let inactive = !element.classList.contains("active");
        setClass(group_selectors, "active", "remove");

        if (inactive) {
            element.classList.toggle("active");
        }
    }


    let activeHm = get_active_element("hm-selector");

    let possibleGenes = combinations[activeHm.value][element.value];

    let dropdown = document.getElementById("select-gene");
    let l, k = dropdown.options.length - 1;

    for (l = k; l >= 0; l--) {
        dropdown.remove(l);
    }

    let m;
    for (m = 0; m < possibleGenes.length; m++) {
        let option = document.createElement("option");
        option.value = possibleGenes[m];
        option.textContent = (m + 1) + ". " + possibleGenes[m];
        dropdown.appendChild(option);
    }

    update_image();
}

function moveGene(delta) {
    let dropdown = document.getElementById("select-gene");
    let options = dropdown.children;
    let i;

    for (i = 0; i < options.length; i++) {
        if (options[i].value === document.getElementById("select-gene").value) {
            break;
        }
    }

    i += delta;

    if (i >= options.length) {
        i = 0;
    }
    if (i < 0) {
        i = options.length - 1;
    }

    dropdown.value = options[i].value;
    update_image();
}

function update_image() {
    let activeHm = get_active_element("hm-selector").value;
    let activeGroup = get_active_element("group-selector").value;
    let selectedGeneName = document.getElementById("select-gene").value;

    let image = document.getElementById("validation-plot");
    image.src = activeHm + "/" + activeGroup + "/" + selectedGeneName + ".png";

    let modal_image = document.getElementById("validation-plot-modal-image");
    modal_image.src = activeHm + "/" + activeGroup + "/" + selectedGeneName + ".png";

    let possibleGenes = combinations[activeHm][activeGroup];
    let modal_caption = document.getElementById("validation-plot-modal-caption");
    modal_caption.innerHTML = (possibleGenes.indexOf(selectedGeneName) + 1) + ". " + selectedGeneName;
}

function get_active_element(className) {
    let elements = document.getElementsByClassName(className);
    let j;

    for (j = 0; j < elements.length; j++) {
        if (elements[j].classList.contains("active")) {
            return elements[j];
        }
    }
}