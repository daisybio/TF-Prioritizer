import {Component, EventEmitter, Input, OnInit, Output} from '@angular/core';

@Component({
  selector: 'accordion-header',
  templateUrl: './header.component.html',
  styleUrls: ['./header.component.css']
})
export class HeaderComponent implements OnInit {

  @Input()
    // @ts-ignore
  title: string;

  @Input()
    // @ts-ignore
  sub: boolean;

  @Output()
  clickEmitter = new EventEmitter();

  constructor() { }

  ngOnInit(): void {
  }

  onClick()
  {
    this.clickEmitter.emit();
  }

  getClasses() {
    return {
      "sub": this.sub
    }
  }
}
