function toggleAccordion(id) {
    let panel = document.getElementById(id);
    if (panel.style.display === "grid") {
        panel.style.display = "none";

        let info = document.getElementById(id + "-info");
        if (info != null) {
            info.style.display = "none";
        }
    } else {
        panel.style.display = "grid";
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

function init_selection(id, combinations) {
    let firstValue;

    if (Array.isArray(combinations)) {
        firstValue = combinations[0];
    } else {
        firstValue = Object.keys(combinations)[0];
    }

    let firstButton = document.getElementById(id + "-0-" + firstValue);
    firstButton.click();
}

function update_selection(source_element, id, combinations) {
    let level = parseInt(source_element.id.split("-")[1]);
    let depth = get_depth(combinations);

    let parentElements = [];

    { // Find activated parent nodes
        for (let i = 0; i < level; i++) {
            let nodes = document.querySelectorAll('[id^="' + id + '-' + i + '"]');

            for (let j = 0; j < nodes.length; j++) {
                if (nodes[j].classList.contains("active")) {
                    parentElements.push(nodes[j].value);
                    break;
                }
            }
        }
    }

    let availableCombinations = combinations;

    for (let i = 0; i < parentElements.length; i++) {
        availableCombinations = availableCombinations[parentElements[i]];
    }

    { // Activate clicked element and deactivate siblings
        let siblings = document.querySelectorAll('[id^="' + id + '-' + level + '"]');
        let allowedSiblings;

        if (Array.isArray(availableCombinations)) {
            allowedSiblings = availableCombinations;
        } else {
            allowedSiblings = Object.keys(availableCombinations);
        }

        for (let i = 0; i < siblings.length; i++) {
            siblings[i].classList.remove("active");
            siblings[i].disabled = !allowedSiblings.includes(siblings[i].value);
        }
        source_element.classList.add("active");
    }

    availableCombinations = availableCombinations[source_element.value];

    { // Fill dropdown if necessary
        if (level === depth - 1 && document.getElementById(id + "-dropdown") != null) {
            let dropdownContent = document.getElementById(id + "-dropdown-content");
            dropdownContent.innerHTML = "";

            let searchbox = document.createElement("input");
            searchbox.type = "text";
            searchbox.placeholder = "Search...";
            searchbox.id = id + "-dropdown-search";
            searchbox.style.width = "100%";
            searchbox.addEventListener("keyup", function () {
                filter_dropdown(id)
            });
            searchbox.classList.add("selector")

            dropdownContent.appendChild(searchbox);

            for (let i = 0; i < availableCombinations.length; i++) {
                let option = document.createElement("button");
                let value = availableCombinations[i];
                option.value = value;
                option.id = id + "-" + 2 + "-" + value;
                option.classList.add("selector");

                let text = value.replace(/\.[^/.]+$/, "");
                if (/^[0-9]+_*/.test(text)) {
                    text = text.substring(text.split("_")[0].length + 1);
                }

                option.textContent = (i + 1) + ". " + text;
                option.addEventListener("click", function () {
                    update_selection(option, id, combinations);
                });

                dropdownContent.appendChild(option);
            }

            let dropdownButton = document.getElementById(id + "-dropdown");
            let dropdownLeftArrow = document.getElementById(id + "-next-option");
            let dropdownRightArrow = document.getElementById(id + "-previous-option");
            let modalLeftArrow = document.getElementById(id + "-modal-leftarrow");
            let modalRightArrow = document.getElementById(id + "-modal-rightarrow");

            dropdownButton.disabled = dropdownLeftArrow.disabled = dropdownRightArrow.disabled = modalLeftArrow.disabled = modalRightArrow.disabled = availableCombinations.length <= 1;
        }
    }

    { // Click first possible child
        if (level < depth) {
            let firstPossibleChildValue;
            if (Array.isArray(availableCombinations)) {
                firstPossibleChildValue = availableCombinations[0];
            } else {
                firstPossibleChildValue = Object.keys(availableCombinations)[0];
            }
            let term = id + "-" + (level + 1) + "-" + firstPossibleChildValue;
            let firstPossibleChildNode = document.getElementById(term);
            firstPossibleChildNode.disabled = false;
            firstPossibleChildNode.click();
        }
    }

    { // Update image if necessary
        if (level === depth) {
            let image = document.getElementById(id + "-image");
            let modal_image = document.getElementById(id + "-modal-image");

            let source = id + "/";

            for (let i = 0; i < parentElements.length; i++) {
                source += parentElements[i] + "/";
            }
            source += source_element.value;

            modal_image.src = image.src = source;

            let modal_caption = document.getElementById(id + "-modal-caption");
            modal_caption.textContent = source_element.textContent;

            let dropdown = document.getElementById(id + "-dropdown");
            if (dropdown != null) {
                dropdown.textContent = source_element.textContent;
                dropdown.value = source_element.value;
            }
        }
    }
}

function move_lowest_level(id, delta, combinations) {
    let options = document.querySelectorAll('[id^="' + id + '-' + get_depth(combinations) + '"]');
    let i;

    for (i = 0; i < options.length; i++) {
        if (options[i].classList.contains("active")) {
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

    options[i].click();

}

function openPanel(id) {
    let panel = document.getElementById(id);
    panel.style.display = "grid";
}

function toggleDropdown(selection) {
    document.getElementById(selection + "-dropdown-content").classList.toggle("show");
}

function filter_dropdown(selection) {
    let searchbox = document.getElementById(selection + "-dropdown-search");
    let term = searchbox.value.toUpperCase();
    let options = document.getElementById(selection + "-dropdown-content").childNodes;

    for (let i = 0; i < options.length; i++) {
        if (options[i] === searchbox) {
            continue;
        }
        if (options[i].textContent.toUpperCase().indexOf(term) > -1) {
            options[i].style.display = "";
        } else {
            options[i].style.display = "none";
        }
    }
}

function add_dropdown_closing(id) {
    document.addEventListener("click", function (e) {
        let content = document.getElementById(id + "-dropdown-content");
        if (!content.contains(e.target) && e.target !== document.getElementById(id + "-dropdown")) {
            content.classList.remove("show");
        }
    });
}

function get_depth(combinations) {
    let i = 0;

    while (!Array.isArray(combinations)) {
        combinations = combinations[Object.keys(combinations)[0]];
        i++;
    }
    return i;
}

function filterTfsByName(term) {
    let tfs = document.getElementsByClassName("tf");

    for (let i = 0; i < tfs.length; i++) {
        if (tfs[i].id.toUpperCase().indexOf(term.toUpperCase()) > -1) {
            tfs[i].style.display = "";
        } else {
            tfs[i].style.display = "none";
        }
    }
}

function filterTfsByTargetGene(term) {
    let tfs = document.getElementsByClassName("tf");

    for (let i = 0; i < tfs.length; i++) {
        let varName = "targetGenes" + tfs[i].id;
        let targetGenes = window[varName];
        let found = false;

        for (let j = 0; j < targetGenes.length; j++) {
            if (targetGenes[j].toUpperCase().indexOf(term.toUpperCase()) > -1) {
                found = true;
                break;
            }
        }

        if (found) {
            tfs[i].style.display = "";
        } else {
            tfs[i].style.display = "none";
        }
    }
}