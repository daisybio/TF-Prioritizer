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

function init_selection(selection) {
    let group_selectors = document.getElementsByClassName(selection + " group-selector");
    select_group(selection, group_selectors[0]);

    select_first_possible_subgroup(selection);

    update_image(selection);
}

function select_first_possible_subgroup(selection) {
    let subgroups = document.getElementsByClassName(selection + " subgroup-selector");
    let i;

    for (i = 0; i < subgroups.length; i++) {
        if (!subgroups[i].disabled) {
            select_subgroup(selection, subgroups[i]);
            break;
        }
    }
}

function select_group(selection, element) {
    let groups = document.getElementsByClassName(selection + " group-selector");
    let subgroups = document.getElementsByClassName(selection + " subgroup-selector")

    if (element.classList.contains("active")) return;

    let inactive = !element.classList.contains("active");
    setClass(groups, "active", "remove");

    if (inactive) {
        element.classList.toggle("active");
    }

    let possible_subgroups = Object.keys(combinations[element.value]);

    let j;

    for (j = 0; j < subgroups.length; j++) {
        subgroups[j].disabled = !possible_subgroups.includes(subgroups[j].value);
    }

    let active_subgroup = get_active_element(selection + " group-selector");

    if (!(typeof active_subgroup === 'undefined')) {
        if (!possible_subgroups.includes(active_subgroup.value)) {
            select_first_possible_subgroup(selection);
        } else {
            select_subgroup(selection, active_subgroup);
        }
    }
}

function select_subgroup(selection, element) {
    let group_selectors = document.getElementsByClassName(selection + " subgroup-selector");

    if (!element.classList.contains("active")) {
        let inactive = !element.classList.contains("active");
        setClass(group_selectors, "active", "remove");

        if (inactive) {
            element.classList.toggle("active");
        }
    }


    let active_group = get_active_element(selection + " group-selector");

    let possible_dropdown_values = combinations[active_group.value][element.value];

    let dropdown = document.getElementById(selection + "-dropdown");
    let l, k = dropdown.options.length - 1;

    for (l = k; l >= 0; l--) {
        dropdown.remove(l);
    }

    let m;
    for (m = 0; m < possible_dropdown_values.length; m++) {
        let option = document.createElement("option");
        option.value = possible_dropdown_values[m];
        option.textContent = (m + 1) + ". " + possible_dropdown_values[m];
        dropdown.appendChild(option);
    }

    update_image(selection);
}

function move_dropdown(selection, delta) {
    let dropdown = document.getElementById(selection + "-dropdown");
    let options = dropdown.children;
    let i;

    for (i = 0; i < options.length; i++) {
        if (options[i].value === dropdown.value) {
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
    update_image(selection);
}

function update_image(selection) {
    let active_group = get_active_element(selection + " group-selector").value;
    let active_subgroup = get_active_element(selection + " subgroup-selector").value;
    let active_dropdown = document.getElementById(selection + "-dropdown").value;

    let image = document.getElementById(selection + "-image");
    image.src = active_group + "/" + active_subgroup + "/" + active_dropdown + ".png";

    let modal_image = document.getElementById(selection + "-modal-image");
    modal_image.src = active_group + "/" + active_subgroup + "/" + active_dropdown + ".png";

    let possibleGenes = combinations[active_group][active_subgroup];
    let modal_caption = document.getElementById(selection + "-modal-caption");
    modal_caption.innerHTML = (possibleGenes.indexOf(active_dropdown) + 1) + ". " + active_dropdown;
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