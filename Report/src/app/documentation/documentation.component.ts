import {Component, OnInit} from '@angular/core';

import {loremIpsum} from "../../assets/information";

@Component({
  selector: 'app-documentation',
  templateUrl: './documentation.component.html',
  styleUrls: ['./documentation.component.css']
})
export class DocumentationComponent implements OnInit {
  loremIpsum: string = loremIpsum;

  constructor() {
  }

  ngOnInit(): void {
  }

}
